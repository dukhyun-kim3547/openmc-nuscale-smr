# =============================================================================
# NuScale VOYGR vs 대형 PWR 핀 셀 비교 분석
# 목적: 설계 파라미터 차이가 중성자 물리에 미치는 영향 정량화
# 비교 항목: 농축도, 냉각재 밀도, 온도 효과
# =============================================================================

import openmc
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# =============================================================================
# 비교할 설계 파라미터
# =============================================================================

designs = {
    'Large PWR\n(Westinghouse AP1000)': {
        'enrichment': 3.1,
        'water_density': 0.71,   # g/cc, 325°C
        'pitch': 1.26,
        'color': 'steelblue',
        'marker': 'o'
    },
    'NuScale VOYGR\n(SMR)': {
        'enrichment': 4.95,
        'water_density': 0.74,   # g/cc, 300°C (더 낮은 온도)
        'pitch': 1.26,
        'color': 'darkorange',
        'marker': 's'
    }
}

# 농축도 범위 비교 (1% ~ 20%)
enrichments = [1.0, 2.0, 3.0, 3.1, 4.0, 4.95, 5.0, 7.0, 10.0, 15.0, 19.75]

def run_pin_cell(enrichment, water_density, pitch, label=""):
    """핀 셀 k-infinity 계산 함수"""

    # 재료 정의
    uo2 = openmc.Material(name=f'UO2_{enrichment}pct')
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

    materials = openmc.Materials([uo2, zircaloy, water])
    materials.export_to_xml()

    # 기하학
    fuel_or = openmc.ZCylinder(r=0.4096)
    clad_or = openmc.ZCylinder(r=0.4750)
    left   = openmc.XPlane(x0=-pitch/2, boundary_type='reflective')
    right  = openmc.XPlane(x0=+pitch/2, boundary_type='reflective')
    bottom = openmc.YPlane(y0=-pitch/2, boundary_type='reflective')
    top    = openmc.YPlane(y0=+pitch/2, boundary_type='reflective')

    fuel_cell  = openmc.Cell(fill=uo2,      region=-fuel_or)
    clad_cell  = openmc.Cell(fill=zircaloy, region=+fuel_or & -clad_or)
    water_cell = openmc.Cell(fill=water,    region=+clad_or & +left & -right & +bottom & -top)

    universe = openmc.Universe(cells=[fuel_cell, clad_cell, water_cell])
    geometry = openmc.Geometry(universe)
    geometry.export_to_xml()

    # 설정
    settings = openmc.Settings()
    settings.run_mode  = 'eigenvalue'
    settings.batches   = 150
    settings.inactive  = 50
    settings.particles = 5000

    bounds = [-0.4096, -0.4096, -1, 0.4096, 0.4096, 1]
    settings.source = openmc.IndependentSource(
        space=openmc.stats.Box(bounds[:3], bounds[3:])
    )
    settings.export_to_xml()

    openmc.run(output=False)

    sp = openmc.StatePoint('statepoint.150.h5')
    keff = sp.keff
    sp.close()

    return keff.nominal_value, keff.std_dev

# =============================================================================
# PART 1: 농축도별 k-infinity 비교 (PWR vs NuScale 냉각재 밀도)
# =============================================================================

print("=" * 55)
print("  NuScale vs PWR: 농축도별 k-infinity 비교")
print("=" * 55)

results = {}
for name, params in designs.items():
    print(f"\n[{name.replace(chr(10), ' ')}] 계산 중...")
    keffs, kerrs = [], []
    for enr in enrichments:
        k, err = run_pin_cell(enr, params['water_density'], params['pitch'])
        keffs.append(k)
        kerrs.append(err)
        print(f"  농축도 {enr:5.2f}% → k = {k:.5f} ± {err:.5f}")
    results[name] = {'keffs': keffs, 'kerrs': kerrs}

# =============================================================================
# PART 2: 주요 운전점 직접 비교
# =============================================================================

print("\n" + "=" * 55)
print("  주요 운전점 비교")
print("=" * 55)
print(f"{'설계':>25} | {'농축도':>6} | {'k-infinity':>10} | {'냉각재밀도':>10}")
print("-" * 60)

key_points = []
for name, params in designs.items():
    enr = params['enrichment']
    k, err = run_pin_cell(enr, params['water_density'], params['pitch'])
    key_points.append((name, enr, k, err, params['water_density']))
    short_name = name.replace('\n', ' ')
    print(f"{short_name:>25} | {enr:>6.2f}% | {k:>8.5f}±{err:.5f} | {params['water_density']:>8.3f} g/cc")

# =============================================================================
# PART 3: 그래프
# =============================================================================

fig, axes = plt.subplots(1, 2, figsize=(13, 5))

# 왼쪽: 농축도별 k 비교 곡선
for name, params in designs.items():
    keffs = results[name]['keffs']
    kerrs = results[name]['kerrs']
    axes[0].errorbar(enrichments, keffs, yerr=kerrs,
                     fmt=params['marker']+'-', color=params['color'],
                     capsize=3, linewidth=2, markersize=6,
                     label=name.replace('\n', ' '))

axes[0].axhline(y=1.0, color='red', linestyle='--', linewidth=1.5, label='Critical (k=1.0)')
axes[0].axvline(x=5.0, color='gray', linestyle=':', linewidth=1.2, label='LEU limit (5%)')

# 각 설계 운전점 표시
for name, params in designs.items():
    enr = params['enrichment']
    idx = enrichments.index(enr) if enr in enrichments else None
    if idx is not None:
        axes[0].scatter([enr], [results[name]['keffs'][idx]],
                       s=150, color=params['color'], zorder=5,
                       edgecolors='black', linewidth=1.5)

axes[0].set_xlabel('U-235 Enrichment (wt%)', fontsize=12)
axes[0].set_ylabel('k-infinity', fontsize=12)
axes[0].set_title('k-infinity vs Enrichment\nPWR vs NuScale VOYGR', fontsize=11)
axes[0].legend(fontsize=9)
axes[0].grid(True, alpha=0.3)

# 오른쪽: 막대 그래프 직접 비교
labels = ['Large PWR\n(3.1%, 0.71 g/cc)', 'NuScale VOYGR\n(4.95%, 0.74 g/cc)']
k_vals = [key_points[0][2], key_points[1][2]]
k_errs = [key_points[0][3], key_points[1][3]]
colors = ['steelblue', 'darkorange']

bars = axes[1].bar(labels, k_vals, yerr=k_errs, color=colors,
                   capsize=5, edgecolor='black', linewidth=0.8,
                   alpha=0.85, width=0.5)
axes[1].axhline(y=1.0, color='red', linestyle='--', linewidth=1.5, label='Critical')
axes[1].set_ylabel('k-infinity', fontsize=12)
axes[1].set_title('Operating Point Comparison\nPWR vs NuScale', fontsize=11)
axes[1].set_ylim(1.2, 1.6)
axes[1].legend(fontsize=10)
axes[1].grid(True, alpha=0.3, axis='y')

for bar, k, err in zip(bars, k_vals, k_errs):
    axes[1].text(bar.get_x() + bar.get_width()/2, k + err + 0.005,
                f'k = {k:.4f}', ha='center', va='bottom', fontsize=10, fontweight='bold')

plt.suptitle('NuScale VOYGR vs Large PWR: Neutron Physics Comparison', fontsize=12)
plt.tight_layout()
plt.savefig('nuscale_vs_pwr.png', dpi=150, bbox_inches='tight')
print("\n✅ 그래프 저장됨: nuscale_vs_pwr.png")
