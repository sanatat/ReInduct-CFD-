# ReInduct-CFD Pro v8 — راهنمای فارسی و ممیزی اجرایی

**تاریخ:** ۵ نوامبر ۲۰۲۶. این بسته حاصل اصلاح v7 است، نه ادعای ساخت یک Digital Twin صنعتی.

## اصلاحات به‌صورت کد

1. **اصلاح جرم مولی محصول:** در `v7/validation/interfacial_law.py` و `v7/openfoam/solver/reactiveFixedBedFoamV7.C` سرعت رشد لایه Ni از چگالی مولی **Ni فلزی** (تقریباً 151,772 mol/m³) به نرخ مصرف NiO تبدیل می‌شود؛ قبلاً از چگالی مولی NiO (تقریباً 87,693 mol/m³) استفاده شده بود. نسبت 1.731 است. تخلخل لایه محصول و سطح واکنش همچنان پارامتر نامعلوم‌اند.
2. **اصلاح دمای گونه‌های گازی:** کسر مولی H₂ با `Tg` نه `Ts` محاسبه می‌شود. ماژول `v8/physics.py` خواص تقریبی مخلوط ایده‌آل H₂/H₂O/N₂ و واحدها را آزمون می‌کند. حلگر OpenFOAM همچنان تقریب چگالی ثابت دارد و **حلگر تراکم‌پذیر کامل نیست**.
3. **تابش دیواره:** sink تابشی خاکستری به مدل C++ اضافه شده؛ استفاده از wall area/volume مدل **تقریبی lumped** است. `v8/thermal_em.py` انتقال حرارت در مدل محورمتقارن را با هدایت finite volume، انتقال گاز–جامد و تابش حل می‌کند.
4. **کوپل دوطرفه میدان و حرارت:** `v8/thermal_em.py` در هر گام تکرار، دمای جامد را به حل الکترومغناطیسی منتقل و توان ژولی را مجدداً محاسبه می‌کند. **این کوپل آفلاین Python است، نه OpenFOAM runtime**.
5. **نگاشت خروجی OpenFOAM:** `v8/openfoam/feedback.py` دما را از فایل ASCII scalar OpenFOAM می‌خواند، از مختصات واقعی سلول‌ها استفاده می‌کند و qEM را با حفظ توان کل بازمی‌نویسد. این رابط با داده‌های مصنوعی تست شده، نه با خروجی واقعی OpenFOAM.
6. **آزمون اصلاح فشار–سرعت:** در فایل C++، تصحیح شار `phiHbyA - pEqn.flux()` و گزارش `maxDivergencePhi` اضافه شده؛ **کامپایل/اجرای واقعی این حلگر هنوز انجام نشده**.
7. **LiCoO₂ واقعی‌تر:** `v8/lco_kinetics.py` و `v8/lco_holdout.py` از مکانیسم diffusion/GB و nucleation/A3 و Eaهای 121 و 95 kJ/mol استفاده می‌کنند. ضرایب سرعت و زمان گذار همچنان مفروض‌اند؛ برازش روی دو دمای 800 و 1000°C و آزمون روی 700 و 1100°C با **چهار نقطه تقریبی همان مقاله** انجام شده، نه آزمایش مستقل.
8. **آزمون مستقل C++17:** `v8/cpp_kinetics_check.cpp` کامپایل و اجرا می‌شود و حفظ موازنه واکنش را با چگالی تصحیح‌شده آزمایش می‌کند.

## اجرای محلی (Linux Mint / Ubuntu)

```bash
sudo apt install python3 python3-venv python3-pip g++ make
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
make v8all
```

شامل ۷۶ آزمون Python، اجرای مدل دوطرفه EM–حرارت، نمودارها و تست مستقل C++.

خروجی: `v8/results/` به‌خصوص `EVIDENCE_MANIFEST.json`, `em_thermal/`, `lco_holdout.json`, `all_pytests.log`, `openfoam_status.txt`.

برای تلاش واقعی OpenFOAM پس از نصب و فعال‌کردن v2606:

```bash
make v8openfoam
```

این فرمان **در محیط حاضر اجرا نشده**، چون `wmake` موجود نیست. می‌توان CI را از GitHub Actions در `.github/workflows/v8-validation.yml` اجرا و لاگ خطاها را اصلاح کرد؛ تعریف workflow به معنی اجرای موفق آن نیست.

## مرز علمی

- ۷۶ تست، **verification محلی** هستند، نه اعتبارسنجی آزمایشگاهی فرایند.
- مدل full plasma/MHD، electron-energy، CFD–DEM و حلگر تراکم‌پذیر چندگونه‌ای هنوز کامل نشده است.
- داده تجربی مستقل راکتور یا منحنی تبدیل واقعی black mass وجود ندارد.
- محیط پروژه برای آزمایش با سیستم لینوکس معمولی طراحی شده، ولی اجرای CFD کامل به نصب نرم‌افزار نیاز دارد.

## منابع علمی

- Hidayat و همکاران، ۲۰۰۹: https://doi.org/10.1007/s11663-008-9212-0
- LiCoO₂/C kinetics، ۲۰۲۳: https://pmc.ncbi.nlm.nih.gov/articles/PMC9821372/
- OpenFOAM 2606 Setup GitHub Action: https://github.com/gerlero/setup-openfoam
