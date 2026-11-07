# ReInduct-CFD Pro v7 | نسخه ۷ برای TU-BRIDGE

**خروجی واقعی:** ۵۸ آزمون Python موفق، تست مستقل C++17 موفق، نگاشت محافظه‌کارانه گرمایش القایی، بررسی همگرایی مکانی، تحلیل holdout روی داده مقاله، موازنه استوکیومتری LiCoO2 در دو ناحیه دمایی.

**خروجی هنوز تأییدنشده:** کامپایل و اجرای OpenFOAM.com v2606. محیط این اجرا `wmake`/`blockMesh` نداشت و دانلود بسته سیستم از طریق DNS ممکن نبود. کد و گردش‌کار GitHub Actions برای اجرای واقعی آماده شده؛ به‌محض اجرای Action باید بر اساس لاگ‌ها اشکالات احتمالی اصلاح شوند. هیچ فایل `VERIFIED.txt` برای OpenFOAM تولید نشده است.

## ۱ — اجرای OpenFOAM
- حلگر جدید: `v7/openfoam/solver/reactiveFixedBedFoamV7.C`
- ابزار استخراج هندسه واقعی سلول‌ها: `v7/openfoam/cellGeometry/exportCellGeometry.C`
- کیس محورمتقارن annular wedge: `v7/openfoam/case`
- دستور یکپارچه: `bash v7/openfoam/Allrun`
- گردش‌کار اجرای واقعی: `.github/workflows/v7-validation.yml`
- وضعیت: نوشته‌شده و از لحاظ ساختاری کنترل شده، اما **compile/run واقعی انجام نشده است**. خطاهای API/مش احتمالی باید از روی لاگ CI برطرف شوند.

## ۲ — معادلات، واحدها و کوپلینگ
در مدل OpenFOAM جدید، `H2` و `H2O` غلظت مولی گاز در منافذ [mol/m³ gas] هستند؛ نه کسر جرمی یا غلظت بی‌بعد. `X` کسر تبدیل جامد و `R` [mol/m³ bulk/s] است. برای تخلخل ثابت ε:

- ∂(ε cH₂)/∂t + ∇·(u_super cH₂) = ∇·(εD∇cH₂) − R
- ∂(ε cH₂O)/∂t + ∇·(u_super cH₂O) = ∇·(εD∇cH₂O) + R
- ∂X/∂t = R/N_NiO ؛ N_NiO [mol/m³ bulk]
- معادله انرژی جامد شامل `qEM − ΔH R` است؛ `qEM` با واحد [W/m³ bulk].

در هر گام زمانی چند تکرار Picard برای سرعت، غلظت و دما تعریف شده است؛ حلقه، همگرایی سراسری کامل را تضمین نمی‌کند. در آزمون مستقل صفر‌بعدی، موازنه جرم و انرژی با خطاهای شناور نزدیک صفر بررسی شده است. این **verification** است، نه validation صنعتی.

## ۳ — نگاشت محافظه‌کارانه میدان القایی
تابع `v7/mapping.py` مدل axisymmetric v6 را حل می‌کند و با مختصات و حجم واقعی سلول‌های صادرشده از OpenFOAM نگاشت می‌کند. قبل از نگاشت، توان مرجع در دامنه استوانه‌ای annular با انتگرال `q(r,z) r dr dtheta dz` محاسبه می‌شود. در انتها منبع سلولی طوری نرمال می‌شود که مجموع `q_i V_i` توان مرجع را حفظ کند.

آزمون آفلاین فعلی روی شبکه ساخته‌شده به‌صورت ریاضی (نه مش خروجی OpenFOAM) اجرا شده است:

| شبکه radial×axial | خطای L2 موضعی پس از نرمال‌سازی | اختلاف توان پس از نرمال‌سازی |
|---|---:|---:|
| ۶×۲۴ | ۰٫۷۴۱٪ | در حد ماشین |
| ۱۲×۴۸ | ۰٫۱۳۸٪ | در حد ماشین |
| ۲۴×۹۶ | ۰٫۱۱۱٪ | در حد ماشین |
| ۳۶×۱۴۴ | ۰٫۰۴۲۸٪ | در حد ماشین |

اعداد از `v7/results/mapping_convergence/mapping_convergence.json` گرفته شده‌اند. محافظه‌کاربودن توان به معنی دقت محلی خودکار نیست؛ به همین علت خطای L2 جداگانه گزارش شده است.

## ۴ — داده علمی مستقل از برازش (در همان مطالعه)
داده ۱۹نقطه‌ای جدول III مقاله Hidayat و همکاران (DOI: 10.1007/s11663-008-9212-0) در `v7/data/nio_hidayat_2009_table_III.csv` قرار دارد. دمای ۸۰۰ درجه کاملاً از آموزش کنار گذاشته شد:

- Arrhenius تک‌رژیمی: MAE = ۰٫۶۳۴ µm/s
- مدل تجربی PCHIP: MAE = ۰٫۵۷۶ µm/s

این ضعف پیش‌بینی در گذار غیرتک‌نواخت دما را آشکار می‌کند. همچنین غلظت H₂ برابر ۲۶٪ به‌صورت جداگانه از آموزش حذف و ارزیابی شد. **این آزمون خارج از نمونه در همان مقاله است**، نه اعتبارسنجی آزمایشگاه مستقل.

قانون جدید اتصال سطح به حجم با فرض سطح مؤثر `a_eff` تعریف شده:

`R_bulk = v_Ni_layer(T,yH2) × (rho_NiO/M_NiO) × a_eff × (1−X)`

که در آن `v` [m/s] و `a_eff` [m²/m³ bulk] است. مقدار a_eff = ۲ فرضی است و با داده راکتور کالیبره نشده. فایل `constant/kineticTable` برای حلگر صادر می‌شود؛ پس مدل NiO فعلاً «تعمیم فرضی نرخ سطحی به نرخ حجمی» است.

برای LiCoO2/C دو وضعیت گزارش تجربی بازنمایی شده‌اند: محصولات اصلی `Li₂CO₃ + Co` در دماهای ۷۰۰–۹۰۰°C و `Li₂O + Co` در ۱۰۰۰–۱۱۰۰°C. واکنش‌های نمونه از لحاظ اتمی موازنه شده‌اند ولی شیمی CO/CO₂ و درصد تبدیل هنوز **فرضی** هستند. هیچ بازده آزمایشگاهی بازیافت باتری از این بخش نتیجه نمی‌شود.

## اجرای کل آزمون‌ها روی Linux Mint

```bash
sudo apt install g++ make python3 python3-venv
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
bash v7/run_v7.sh
```

برای **اجرای واقعی OpenFOAM 2606** باید نسخه سازگار نصب و `wmake`, `blockMesh`, `checkMesh` در مسیر شِل باشند:

```bash
bash v7/openfoam/Allrun
```

یا همین بسته را در repository خصوصی GitHub بارگذاری و workflow `v7-validation.yml` را از Actions اجرا کن. فقط پس از موفقیت Job مربوط به OpenFOAM می‌توان ادعای کامپایل و اجرا داشت.

## محدوده ادعاهای قابل دفاع در CV
**مجاز:** numerical conservation verification, experimental kinetic holdout analysis (NiO), C++17 scientific calculations, axisymmetric EM heat-source mapping with grid-refinement, conservative reactive transport source formulation, literature-informed battery mass balances.

**فعلاً غیرمجاز:** successfully executed OpenFOAM reactive CFD, validated industrial battery reactor digital twin, experimentally calibrated NiO packed-bed kinetics, full plasma–MHD coupling.
