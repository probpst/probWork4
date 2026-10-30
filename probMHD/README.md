# probMHD — full-MHD solver for tokamak plasmas (probCFD_부산 family)

**probMHD**는 probCFD_부산의 코드 양식(그룹형 name list, PLOT3D 입출력, MPI, 정지 기준과
`Reason_for_stop.md`)을 따르는 **완전(Full) 압축성 점성·저항 단일유체 MHD 솔버**입니다.
구조 곡선좌표 격자 위의 유한체적법으로, 토카막(KSTAR 규모) 내부 플라즈마의 거시 불안정성을
계산하기 위해 만들어졌습니다. probCFD_부산의 **2차 ENO(minmod) 대류 스킴**(`convection.F90`의
`upw`/`minmod`)을 그대로 이식했고, 난류 모델 **스위치**(`SW_Turb_MHD`)로 모델을 끄고 켤 수 있습니다.

| 항목 | 내용 |
|---|---|
| 방정식 | 압축성 저항·점성 MHD: 연속, 운동량, 에너지(전체 에너지 또는 내부 에너지), 유도 방정식 |
| 이산화 | 셀 중심 유한체적(곡선좌표, 단일 블록), 면 자속의 **제약 수송(CT)**으로 ∇·B = 0 유지 |
| 대류 | 리만 솔버 Rusanov / HLL / **HLLD**; 재구성 1차 upwind / 2차 upwind / QUICK / **2차 ENO(probCFD)** / 3차 ENO |
| 시간 | 명시적 Euler / SSP-RK2 / **SSP-RK3**, CFL·VN 자동 시간 간격 |
| 비이상 항 | 저항(Spitzer 선택), 점성, **비등방 열전도**(χ∥, χ⊥), 입자 확산, 옴 가열, 가열·제동복사 소스 |
| 난류 스위치 | `SW_Turb_MHD` 0 끔 / 1 MHD Smagorinsky(ν_sgs, η_sgs) / 3 이상 수송(χ_anom, Bohm) / 2 동적 모델(예약) |
| 평형 | Cerfon–Freidberg Solov'ev 해석 평형(KSTAR 격자 생성기와 동일 구성), 원통 테어링 평형, 표준 시험 |
| 토카막 보조 | 평형 잔차 차감(`SW_Eq_Subtract`), 저β용 내부 에너지 방정식(`SW_Energy`), 토로이달 쐐기 회전 주기성 |
| 격자 | PLOT3D ASCII(probCFD) 또는 Fortran unformatted(KSTAR 파일), 생성 박스; j/k 주기(중복 절점 유무 모두) |
| 병렬 | MPI 3차원 분할(단일 블록), 고스트 교환, rank 0 취합 출력 |

## 폴더 구성

| 경로 | 내용 |
|---|---|
| `src/*.F90`, `Makefile` | 솔버 소스(13개 모듈). `make` (gfortran + OpenMPI), `make debug` |
| `probMHD_template.nml` | **모든 입력 변수**가 주석으로 설명된 name list 템플릿 |
| `cases/grids/` | KSTAR 격자 축소본(33×60 축대칭 쐐기, 99×150, 33×60×20 전체 토러스, 90° 쐐기), 원통 격자 |
| `cases/KSTAR_equilibrium_2D/` | KSTAR 평형 유지 시험(축대칭) |
| `cases/KSTAR_3D_torus/` | 전체 토러스 3-D, m=3/n=1 섭동, MPI 4 rank |
| `cases/KSTAR_turbulence_switch/` | 같은 섭동 평형에서 난류 스위치 0/1/3 비교 |
| `verification/` | Brio–Wu, Orszag–Tang, Alfvén 파 수렴, 테어링 모드 성장률 (스크립트와 결과 그림 포함) |
| `tools/` | PLOT3D 판독(`readp3d.py`), KSTAR 격자 축소/쐐기 생성(`coarsen_plot3d.py`), 원통 격자(`make_cylinder_grid.py`), 모니터 플롯, Solov'ev 계수 교차검증 |

## 빌드와 실행

```bash
cd probMHD && make
cd cases/KSTAR_equilibrium_2D && cp ../../probMHD . && mpirun -np 1 ./probMHD
cd ../KSTAR_3D_torus && cp ../../probMHD . && mpirun -np 4 ./probMHD
python3 ../../tools/plot_monitor.py Output_data/KSTAR_3D_monitor.dat
```

출력: `Volume_solutions_grids/<case>_cells.grd`(셀 중심 격자 [m]), `<case>_NNNNNNNN.p3d`(17개 변수, `.nam`에
이름과 SI 단위), `Output_data/<case>_monitor.dat`(에너지, 최대 속도, ∇·B, 섭동 에너지 이력),
`<case>_slice_NNNNNNNN.dat`(k 중간면 ASCII), `<case>_NNNNNNNN.rsto`(재시작).

## 무차원화와 단위

입력과 출력은 SI입니다(`Units = 1`이면 코드 단위). 코드 단위: 길이 `Char_L`, 자기장 `Char_B`,
밀도 `Char_Rho`; v_A0 = B0/√(μ0 ρ0), t_A0 = L0/v_A0, p0 = B0²/μ0, 저항 η* = η/(μ0 L0 v_A0) = 1/S.
온도는 T* = p*/ρ*이고 T[K] = m_i p0/(2 k_B ρ0)·T* (수소동위원소, Z = 1, T_i = T_e).

## 토카막 계산에서 꼭 알아야 할 세 가지

1. **저β 강성.** 토카막은 p ≪ B²/2 (KSTAR 격자 평형에서 β ~ 10⁻³)이므로 전체 에너지에서 압력을
   빼내면 자기 에너지의 10⁻⁴ 오차가 압력의 10 % 오차가 됩니다. 토카막에는 `SW_Energy = 2`
   (내부 에너지 방정식)를 쓰십시오. 충격파가 있는 표준 시험은 `SW_Energy = 1`(보존형)입니다.
2. **이산 평형 잔차.** 해석 평형도 거친 격자에서는 J×B − ∇p의 이산 잔차가 남아 평형이
   흐릅니다(33×60 격자에서 최대 |v| ~ 0.06 v_A 로 포화). `SW_Eq_Subtract = 1`은 초기 상태의
   이산 잔차를 저장해 빼므로 평형이 정확한 고정점이 되고, 섭동만 발전합니다. 격자를 세분하면
   잔차는 2차로 줄어듭니다(`cases/KSTAR_equilibrium_2D`의 33×60 대 99×150 비교).
3. **토로이달 쐐기 기하.** 데카르트 벡터로 표현한 토러스 쐐기는 현(chord) 오차가 (Δφ)²에 비례합니다.
   축대칭 계산은 `tools/coarsen_plot3d.py --nphi 1 --close --wedge-deg 0.5`로 얇은 쐐기를
   만들고 `Periodic_rot = 0, 0, 0.5`를 쓰십시오. 전체 토러스(KSTAR 3.6° 평면)는 그대로 읽을 수 있습니다.

## 검증 결과 요약

`verification/` 안의 README.md와 그림을 보십시오(Brio–Wu 충격관, Orszag–Tang 와류의 에너지 보존과
∇·B, 원편광 Alfvén 파의 수렴 차수(ENO2 ≈ 1.5, ENO3 ≈ 3.0), 원통 테어링 모드의 S 스케일링,
KSTAR 평형 유지, 난류 스위치 비교).

## 현재 상태에서 확인된 수치적 한계 (중요)

* **격자 규모 수치 저항.** 전단 자기장(전류가 흐르는 평형)에서는 CT 스킴의 수치 저항(~c·Δx)이 격자
  규모(파장 4셀)의 저항성 모드를 키웁니다. 원통 테어링 시험에서 η = 0인데도 γ t_A ≈ 0.3(28×32×8),
  0.65(56×64×8)로 성장했고, 데카르트 전단 슬랩(`IC_type = 7`)에서도 Rusanov 플럭스로 같은 현상이
  납니다(`CT_scheme = 0`이나 균일 자기장에서는 없음). 물리적 S 스케일링은 저항층(δ ~ S^-2/5 a)을
  여러 셀로 분해하는 격자(S = 10³~10⁴에서 반경 방향 수백 셀)에서만 측정할 수 있으며, 이 환경(유효
  1코어)에서는 돌리지 못했습니다. `verification/tearing`의 스크립트는 그 계산을 위해 준비되어 있습니다.
* **자유 평형의 붕괴.** `SW_Eq_Subtract = 0`인 KSTAR 평형은 내부 구멍 경계의 첫 셀 고리에서 자기장과
  밀도를 잃으며 2~3 t_A 안에 무너집니다(99×150에서도 동일). 토카막 계산은 반드시 `SW_Eq_Subtract = 1`과
  `SW_Energy = 2`를 함께 쓰십시오. 섭동을 준 차감 계산(`cases/KSTAR_turbulence_switch`)은 난류 모델의
  유무에 따른 차이를 보여 주지만, 위의 수치 저항 때문에 1e-3 v_A 섭동이 수 t_A 안에 ~0.3 v_A까지 자랍니다.
* 따라서 "점화가 꺼지는 원인"을 이 코드로 판단하려면 (1) 더 촘촘한 격자(또는 3차 ENO)와 (2) 명시적
  저항·점성(또는 난류 모델)으로 수치 저항을 물리 값 아래로 내리는 수렴 확인이 먼저 필요합니다.

## 한계와 다음 단계

* 명시적 시간 전진: Alfvén CFL이 시간 간격을 정합니다(KSTAR 3M 격자에서 Δt ≈ 10⁻¹⁰ s).
  느린 저항 불안정성의 장시간 계산에는 반음해(semi-implicit) 연산자가 필요하며 2단계 과제입니다.
* 단일 격자 블록(MPI 분할). 다중 블록 인터페이스는 probCFD의 b2b 엔진을 이식할 자리로 남겨 두었습니다.
* Solov'ev 평형은 X-점이 없는 코어 격자이며 q = 1, 2 면이 없습니다(KSTAR 격자 README 참조).
  실제 방전의 q 분포는 EFIT 평형에서 격자와 ψ(R,Z)를 만들어 넣어야 합니다.
* 동적 Smagorinsky(`SW_Turb_MHD = 2`)와 Hall 항은 예약된 스위치입니다.
