# =============================================================================
# SMR 소형 노심 임계 크기 탐색 (Critical Size Search)
# 목적: 노심 반지름 변화에 따른 k-effective 계산
#       → 최소 임계 크기(Minimum Critical Size) 결정
# 물리: 노심이 작을수록 중성자 누설 증가 → k-eff 감소
# =============================================================================

import openmc
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# =============================================================================
# NuScale VOYGR 재료 파라미터
# =============================================================================

def make_materials(enrichment=4.95, water_density=0.74):
    uo2 = openmc.Material(name='UO2')
    uo2.add_nuclide('U235', enrichment,           percent_type='wo')
    uo2.add_nuclide('U238', 100-enrichment-13.5,  percent_type='wo')
    uo2.add_element('O',    13.5,                 percent_type='wo')
    uo2.set_density('g/cm3', 10.97)

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

    return uo2, zircaloy, water

def run_finite_core(radius, height, enrichment=4.95, water_density=0.74):
    """
    유한 원통형 노심 k-effective 계산
    radius: 노심 반지름 (cm)
    height: 노심 높이 (cm)
    반사 경계 없음 → 진공 경계 → 중성자 누설 발생
    """

    uo2, zircaloy, water = make_materials(enrichment, water_density)
    materials = openmc.Materials([uo2, zircaloy, water])
    materials.export_to_xml()

    # 핀 셀 유니버스
    pitch = 1.26
    fuel_or = openmc.ZCylinder(r=0.4096)
    clad_or = openmc.ZCylinder(r=0.4750)

    fuel_cell  = openmc.Cell(fill=uo2,      region=-fuel_or)
    clad_cell  = openmc.Cell(fill=zircaloy, region=+fuel_or & -clad_or)
    water_cell = openmc.Cell(fill=water,    region=+clad_or)
    pin_universe = openmc.Universe(cells=[fuel_cell, clad_cell, water_cell])

    # 격자 크기 결정 (반지름에 맞게 핀 수 계산)
    n_pins = max(3, int(2 * radius / pitch))
    if n_pins % 2 == 0:
        n_pins += 1  # 홀수로 맞추기 (중앙 핀 존재)

    half = n_pins * pitch / 2

    # 정사각형 격자로 근사 (원통형 노심 근사)
    lattice = openmc.RectLattice()
    lattice.pitch      = (pitch, pitch)
    lattice.lower_left = (-half, -half)
    lattice.universes  = [[pin_universe]*n_pins for _ in range(n_pins)]

    # 노심 경계: 진공(vacuum) → 중성자 누설 허용
    core_cyl  = openmc.ZCylinder(r=half,      boundary_type='vacuum')
    core_top  = openmc.ZPlane(z0=+height/2,   boundary_type='vacuum')
    core_bot  = openmc.ZPlane(z0=-height/2,   boundary_type='vacuum')

    core_cell = openmc.Cell(fill=lattice,
                            region=-core_cyl & -core_top & +core_bot)
    root_universe = openmc.Universe(cells=[core_cell])
    geometry = openmc.Geometry(root_universe)
    geometry.export_to_xml()

    # 계산 설정
    settings = openmc.Settings()
    settings.run_mode  = 'eigenvalue'
    settings.batches   = 120
    settings.inactive  = 40
    settings.particles = 5000

    bounds = [-half*0.8, -half*0.8, -height/2*0.8,
               half*0.8,  half*0.8,  height/2*0.8]
    settings.source = openmc.IndependentSource(
        space=openmc.stats.Box(bounds[:3], bounds[3:])
    )
    settings.export_to_xml()

    openmc.run(output=False)

    sp = openmc.StatePoint('statepoint.120.h5')
    keff = sp.keff
    leakage = float(sp.global_tallies[2][2])  # leakage fraction
    sp.close()

    return keff.nominal_value, keff.std_dev, n_pins

# =============================================================================
# 노심 크기별 계산
# 높이/직경 비율 = 1.0 (최적 임계 원통 근사)
# =============================================================================

# 반지름 범위: 20cm ~ 120cm
radii = [20, 30, 40, 50, 60, 70, 85, 100, 120]

print("=" * 60)
print("  NuScale SMR 노심 임계 크기 탐색")
print("=" * 60)
print(f"{'반지름(cm)':>10} | {'높이(cm)':>8} | {'핀 수':>6} | {'k-eff':>10} | {'불확도':>8}")
print("-" * 55)

keff_vals = []
keff_errs = []
actual_radii = []

for r in radii:
    h = r * 2.0  # 높이 = 직경 (최적 임계 원통)
    k, err, n = run_finite_core(r, h)
    keff_vals.append(k)
    keff_errs.append(err)
    actual_radii.append(r)
    status = "✅ 임계가능" if k > 1.0 else "❌ 미임계"
    print(f"{r:>10} | {h:>8.0f} | {n:>4}×{n:<2} | {k:>8.5f}±{err:.5f} | {status}")

# =============================================================================
# 최소 임계 반지름 추정 (선형 보간)
# =============================================================================

critical_r = None
for i in range(len(keff_vals)-1):
    if keff_vals[i] < 1.0 and keff_vals[i+1] >= 1.0:
        # 선형 보간
        r1, r2 = actual_radii[i], actual_radii[i+1]
        k1, k2 = keff_vals[i], keff_vals[i+1]
        critical_r = r1 + (1.0 - k1) / (k2 - k1) * (r2 - r1)
        break
    elif keff_vals[i] >= 1.0 and keff_vals[i-1] < 1.0:
        r1, r2 = actual_radii[i-1], actual_radii[i]
        k1, k2 = keff_vals[i-1], keff_vals[i]
        critical_r = r1 + (1.0 - k1) / (k2 - k1) * (r2 - r1)
        break

print(f"\n{'='*55}")
if critical_r:
    print(f"  최소 임계 반지름: ~{critical_r:.1f} cm ({critical_r*2/100:.2f} m 직경)")
    print(f"  최소 임계 높이:   ~{critical_r*2:.1f} cm ({critical_r*4/100:.2f} m)")
print(f"  NuScale 실제 노심: 반지름 ~85cm, 높이 ~200cm")
print(f"{'='*55}")

# =============================================================================
# 그래프
# =============================================================================

plt.figure(figsize=(9, 5))
plt.errorbar(actual_radii, keff_vals, yerr=keff_errs,
             fmt='o-', color='darkorange', capsize=4,
             linewidth=2, markersize=7, label='NuScale VOYGR (4.95% enr.)')
plt.axhline(y=1.0, color='red', linestyle='--',
            linewidth=1.5, label='Critical (k=1.0)')

if critical_r:
    plt.axvline(x=critical_r, color='purple', linestyle=':',
                linewidth=1.5, label=f'Min. critical radius (~{critical_r:.0f} cm)')

plt.axvline(x=85, color='green', linestyle='-.',
            linewidth=1.5, label='NuScale actual (~85 cm)')

plt.fill_between(actual_radii, keff_vals, 1.0,
                 where=[k > 1.0 for k in keff_vals],
                 alpha=0.1, color='green', label='Supercritical region')

plt.xlabel('Core Radius (cm)', fontsize=12)
plt.ylabel('k-effective', fontsize=12)
plt.title('NuScale VOYGR: k-effective vs Core Size\n(Vacuum boundary — neutron leakage included)',
          fontsize=11)
plt.legend(fontsize=9)
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('critical_size.png', dpi=150, bbox_inches='tight')
print("\n✅ 그래프 저장됨: critical_size.png")
