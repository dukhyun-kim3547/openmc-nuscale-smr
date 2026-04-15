# =============================================================================
# HALEU 연소 계산 (HALEU Depletion Analysis)
# 목적: NuScale VOYGR(4.95% HALEU) vs 대형 PWR(3.1%) 연료 수명 비교
# 핵심: 높은 농축도 → 더 긴 연료 수명 → SMR 경제성 근거
# =============================================================================

import openmc
import openmc.deplete
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import math

chain_path = '/home/kevin/miniconda3/envs/openmc/lib/python3.11/site-packages/openmc_data/depletion/chain_endf_b8.0_pwr.xml'

def build_model(enrichment, water_density, label):
    """핀 셀 모델 생성"""

    pitch  = 1.26
    fuel_r = 0.4096

    uo2 = openmc.Material(name=f'UO2_{label}')
    uo2.add_nuclide('U235', enrichment,           percent_type='wo')
    uo2.add_nuclide('U238', 100-enrichment-13.5,  percent_type='wo')
    uo2.add_element('O',    13.5,                 percent_type='wo')
    uo2.set_density('g/cm3', 10.97)
    uo2.depletable = True
    uo2.volume = math.pi * fuel_r**2 * 1.0

    zircaloy = openmc.Material(name='Zircaloy-4')
    zircaloy.add_element('Zr', 98.0, percent_type='wo')
    zircaloy.add_element('Sn',  1.5, percent_type='wo')
    zircaloy.add_element('Fe',  0.2, percent_type='wo')
    zircaloy.add_element('Cr',  0.1, percent_type='wo')
    zircaloy.set_density('g/cm3', 6.56)

    water = openmc.Material(name='Water')
    water.add_nuclide('H1',  2.0, percent_type='ao')
    water.add_nuclide('O16', 1.0, percent_type='ao')
    water.set_density('g/cm3', water_density)
    water.add_s_alpha_beta('c_H_in_H2O')

    materials = openmc.Materials([uo2, zircaloy, water])

    fuel_or = openmc.ZCylinder(r=fuel_r)
    clad_or = openmc.ZCylinder(r=0.4750)
    left   = openmc.XPlane(x0=-pitch/2, boundary_type='reflective')
    right  = openmc.XPlane(x0=+pitch/2, boundary_type='reflective')
    bottom = openmc.YPlane(y0=-pitch/2, boundary_type='reflective')
    top    = openmc.YPlane(y0=+pitch/2, boundary_type='reflective')

    fuel_cell  = openmc.Cell(name='Fuel',      fill=uo2,      region=-fuel_or)
    clad_cell  = openmc.Cell(name='Clad',      fill=zircaloy, region=+fuel_or & -clad_or)
    water_cell = openmc.Cell(name='Moderator', fill=water,    region=+clad_or & +left & -right & +bottom & -top)

    universe = openmc.Universe(cells=[fuel_cell, clad_cell, water_cell])
    geometry = openmc.Geometry(universe)

    settings = openmc.Settings()
    settings.run_mode  = 'eigenvalue'
    settings.batches   = 100
    settings.inactive  = 30
    settings.particles = 5000

    bounds = [-fuel_r, -fuel_r, -1, fuel_r, fuel_r, 1]
    settings.source = openmc.IndependentSource(
        space=openmc.stats.Box(bounds[:3], bounds[3:])
    )

    return openmc.Model(geometry=geometry, materials=materials, settings=settings)

# =============================================================================
# 두 설계 비교
# =============================================================================

designs = {
    'PWR (3.1%)':        {'enrichment': 3.1,  'water_density': 0.71, 'color': 'steelblue'},
    'NuScale HALEU (4.95%)': {'enrichment': 4.95, 'water_density': 0.74, 'color': 'darkorange'},
}

# 연소 스텝 (6스텝, 약 3년)
timesteps_days = [1, 30, 180, 360, 720, 1080]
timesteps_seconds = [d * 86400 for d in timesteps_days]
power = 174.0  # W

all_results = {}

for label, params in designs.items():
    print(f"\n{'='*55}")
    print(f"  {label} 연소 계산 시작")
    print(f"{'='*55}")

    # 결과 파일 이름 구분
    result_file = f"depletion_{label.split('(')[0].strip().replace(' ','_')}.h5"

    model = build_model(params['enrichment'], params['water_density'], label)

    chain = openmc.deplete.Chain.from_xml(chain_path)
    operator = openmc.deplete.CoupledOperator(model, chain)

    integrator = openmc.deplete.PredictorIntegrator(
        operator, timesteps_seconds, power=power, timestep_units='s'
    )
    integrator.integrate()

    # 결과 파일 이름 변경
    import shutil
    shutil.move('depletion_results.h5', result_file)

    # 결과 읽기
    results = openmc.deplete.Results(result_file)
    time_k, k_data = results.get_keff()
    time_days = time_k / 86400

    keff_vals = k_data[:, 0]
    keff_errs = k_data[:, 1]

    _, u235  = results.get_atoms('1', 'U235')
    _, pu239 = results.get_atoms('1', 'Pu239')

    u235_ratio  = u235  / u235[0]  * 100
    pu239_ratio = pu239 / u235[0]  * 100

    all_results[label] = {
        'time_days':   time_days,
        'keff_vals':   keff_vals,
        'keff_errs':   keff_errs,
        'u235_ratio':  u235_ratio,
        'pu239_ratio': pu239_ratio,
        'color':       params['color']
    }

    # 임계 수명 계산
    critical_life = None
    for i in range(len(keff_vals)-1):
        if keff_vals[i] >= 1.0 and keff_vals[i+1] < 1.0:
            t1, t2 = time_days[i], time_days[i+1]
            k1, k2 = keff_vals[i], keff_vals[i+1]
            critical_life = t1 + (1.0-k1)/(k2-k1) * (t2-t1)
            break

    print(f"\n  결과 요약:")
    print(f"  {'시간(일)':>8} | {'k-eff':>8} | {'U235(%)':>8} | {'Pu239(%)':>9}")
    print(f"  {'-'*42}")
    for i in range(len(time_days)):
        print(f"  {time_days[i]:>8.1f} | {keff_vals[i]:>8.5f} | "
              f"{u235_ratio[i]:>8.2f} | {pu239_ratio[i]:>9.4f}")

    if critical_life:
        print(f"\n  → 임계 수명: ~{critical_life:.0f}일 ({critical_life/365:.1f}년)")

# =============================================================================
# 비교 그래프
# =============================================================================

fig, axes = plt.subplots(1, 3, figsize=(15, 5))

for label, res in all_results.items():
    # k-infinity vs 시간
    axes[0].errorbar(res['time_days'], res['keff_vals'], yerr=res['keff_errs'],
                     fmt='o-', color=res['color'], capsize=3,
                     linewidth=2, markersize=5, label=label)

    # U-235 감소
    axes[1].plot(res['time_days'], res['u235_ratio'],
                 'o-', color=res['color'], linewidth=2, markersize=5, label=label)

    # Pu-239 생성
    axes[2].plot(res['time_days'], res['pu239_ratio'],
                 'o-', color=res['color'], linewidth=2, markersize=5, label=label)

axes[0].axhline(y=1.0, color='red', linestyle='--', linewidth=1.5, label='Critical')
axes[0].set_xlabel('Time (days)')
axes[0].set_ylabel('k-infinity')
axes[0].set_title('k-infinity vs Burnup\nPWR vs NuScale HALEU')
axes[0].legend(fontsize=9)
axes[0].grid(True, alpha=0.3)

axes[1].set_xlabel('Time (days)')
axes[1].set_ylabel('U-235 remaining (%)')
axes[1].set_title('U-235 Depletion\nPWR vs NuScale HALEU')
axes[1].legend(fontsize=9)
axes[1].grid(True, alpha=0.3)

axes[2].set_xlabel('Time (days)')
axes[2].set_ylabel('Pu-239 / initial U-235 (%)')
axes[2].set_title('Pu-239 Production\nPWR vs NuScale HALEU')
axes[2].legend(fontsize=9)
axes[2].grid(True, alpha=0.3)

plt.suptitle('HALEU vs PWR Fuel: Depletion Comparison', fontsize=13)
plt.tight_layout()
plt.savefig('haleu_depletion.png', dpi=150, bbox_inches='tight')
print("\n✅ 그래프 저장됨: haleu_depletion.png")
