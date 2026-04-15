# =============================================================================
# 수동 냉각 안전성 분석 (DHRS: Decay Heat Removal System)
# 목적: 냉각재 온도 상승에 따른 반응도 변화 계산
#       → 부의 온도계수(Negative Temperature Coefficient) 확인
#       → NuScale 수동 안전성의 물리적 근거 수치화
# =============================================================================

import openmc
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# IAPWS-IF97 물 상태방정식 (온도 → 밀도 변환)
# 설치: pip install iapws
try:
    from iapws import IAPWS97
    USE_IAPWS = True
    print("✅ IAPWS 사용 가능 — 정확한 물 밀도 계산")
except ImportError:
    USE_IAPWS = False
    print("⚠️  IAPWS 없음 — 선형 근사 사용")

def water_density(temp_C, pressure_MPa=15.5):
    """온도(°C)와 압력(MPa)으로 냉각수 밀도(g/cc) 계산"""
    if USE_IAPWS:
        state = IAPWS97(T=temp_C + 273.15, P=pressure_MPa)
        return state.rho / 1000.0  # kg/m³ → g/cc
    else:
        # 선형 근사 (PWR 운전 범위)
        return 0.74 - 0.002 * (temp_C - 300)

# =============================================================================
# 분석할 온도 범위
# 정상운전: 300°C
# DHRS 작동: 온도 상승 → 최대 ~350°C
# 냉온정지: ~100°C
# =============================================================================

temperatures = [100, 150, 200, 250, 280, 300, 310, 320, 330, 340]

def run_pin_cell(temp_C, pressure_MPa=15.5):
    """온도별 핀 셀 k-infinity 계산"""

    density = water_density(temp_C, pressure_MPa)

    uo2 = openmc.Material(name='UO2')
    uo2.add_nuclide('U235', 4.95,          percent_type='wo')
    uo2.add_nuclide('U238', 100-4.95-13.5, percent_type='wo')
    uo2.add_element('O',    13.5,          percent_type='wo')
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
    water.set_density('g/cm3', density)
    if temp_C < 340:
        water.add_s_alpha_beta('c_H_in_H2O')

    materials = openmc.Materials([uo2, zircaloy, water])
    materials.export_to_xml()

    pitch = 1.26
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

    return keff.nominal_value, keff.std_dev, density

# =============================================================================
# 온도별 계산
# =============================================================================

print("=" * 65)
print("  NuScale DHRS: 온도별 k-infinity 및 온도계수 계산")
print("=" * 65)
print(f"{'온도(°C)':>8} | {'밀도(g/cc)':>10} | {'k-infinity':>12} | {'불확도':>8}")
print("-" * 50)

keff_vals  = []
keff_errs  = []
densities  = []

for T in temperatures:
    k, err, rho = run_pin_cell(T)
    keff_vals.append(k)
    keff_errs.append(err)
    densities.append(rho)
    print(f"{T:>8} | {rho:>10.4f} | {k:>12.5f} | ±{err:.5f}")

# =============================================================================
# 온도계수 계산
# MTC (Moderator Temperature Coefficient) [pcm/°C]
# MTC = (1/k) * dk/dT * 10^5
# 부의 값 → 온도 상승 시 반응도 감소 → 자동 안전
# =============================================================================

print(f"\n{'='*55}")
print("  온도계수 (MTC) 계산")
print(f"{'='*55}")

mtc_values = []
mtc_temps  = []

for i in range(1, len(temperatures)-1):
    dT = temperatures[i+1] - temperatures[i-1]
    dk = keff_vals[i+1] - keff_vals[i-1]
    mtc = (dk / keff_vals[i]) / dT * 1e5  # pcm/°C
    mtc_values.append(mtc)
    mtc_temps.append(temperatures[i])
    print(f"  {temperatures[i]}°C: MTC = {mtc:+.2f} pcm/°C")

avg_mtc = np.mean(mtc_values)
print(f"\n  평균 MTC = {avg_mtc:+.2f} pcm/°C")
if avg_mtc < 0:
    print("  → 부의 온도계수 확인! NuScale 수동 안전성 검증 ✅")

# 정상운전(300°C) 대비 반응도 변화
k_ref = keff_vals[temperatures.index(300)]
print(f"\n  정상운전(300°C) k = {k_ref:.5f} 기준 반응도 변화:")
for T, k in zip(temperatures, keff_vals):
    drho = (k - k_ref) / (k * k_ref) * 1e5  # pcm
    marker = " ← 정상운전" if T == 300 else ""
    print(f"  {T:>5}°C: Δρ = {drho:+8.1f} pcm{marker}")

# =============================================================================
# 그래프
# =============================================================================

fig, axes = plt.subplots(1, 3, figsize=(15, 5))

# k-infinity vs 온도
axes[0].errorbar(temperatures, keff_vals, yerr=keff_errs,
                 fmt='o-', color='steelblue', capsize=3,
                 linewidth=2, markersize=6)
axes[0].axhline(y=1.0, color='red', linestyle='--',
                linewidth=1.5, label='Critical (k=1.0)')
axes[0].axvline(x=300, color='green', linestyle=':',
                linewidth=1.5, label='Normal operation (300°C)')
axes[0].set_xlabel('Coolant Temperature (°C)', fontsize=11)
axes[0].set_ylabel('k-infinity', fontsize=11)
axes[0].set_title('k-infinity vs Temperature\n(NuScale VOYGR)', fontsize=11)
axes[0].legend(fontsize=9)
axes[0].grid(True, alpha=0.3)

# 냉각재 밀도 vs 온도
axes[1].plot(temperatures, densities, 's-', color='darkorange',
             linewidth=2, markersize=6)
axes[1].axvline(x=300, color='green', linestyle=':',
                linewidth=1.5, label='Normal operation')
axes[1].set_xlabel('Coolant Temperature (°C)', fontsize=11)
axes[1].set_ylabel('Coolant Density (g/cc)', fontsize=11)
axes[1].set_title('Coolant Density vs Temperature\n(IAPWS-IF97 or Linear)', fontsize=11)
axes[1].legend(fontsize=9)
axes[1].grid(True, alpha=0.3)

# MTC vs 온도
axes[2].plot(mtc_temps, mtc_values, '^-', color='purple',
             linewidth=2, markersize=6)
axes[2].axhline(y=0, color='red', linestyle='--',
                linewidth=1.5, label='Zero MTC')
axes[2].fill_between(mtc_temps, mtc_values, 0,
                     where=[m < 0 for m in mtc_values],
                     alpha=0.15, color='green', label='Negative (safe)')
axes[2].set_xlabel('Temperature (°C)', fontsize=11)
axes[2].set_ylabel('MTC (pcm/°C)', fontsize=11)
axes[2].set_title('Moderator Temperature Coefficient\n(Negative = Inherently Safe)', fontsize=11)
axes[2].legend(fontsize=9)
axes[2].grid(True, alpha=0.3)

plt.suptitle('NuScale VOYGR DHRS Safety Analysis:\nPassive Cooling via Negative Temperature Coefficient',
             fontsize=12)
plt.tight_layout()
plt.savefig('dhrs_analysis.png', dpi=150, bbox_inches='tight')
print("\n✅ 그래프 저장됨: dhrs_analysis.png")
