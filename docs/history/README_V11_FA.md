# ReInduct-CFD Pro v11

## چه چیزهایی به v10 اضافه شد؟

**۱. راکتور LiCoO₂/C واقعی در کد:** واکنش‌های نماینده و موازنه‌شده، موجودی مستقل فاز جامد و گونه‌های گاز CO/CO2/N2، انتقال گاز و دما، مدل مقاومت Ergun، تابش و انتقال بین‌فازی، حسابداری ورودی/خروجی عناصر. کد: `v11/integrated_reactor.py`.

**۲. کوپل دوطرفه با گرمایش القایی:** میدان الکترومغناطیس محورمتقارن از `v6/axisymmetric_em.py` خوانده و با دمای جامد محاسبه‌شده مجدداً حل می‌شود؛ نگاشت انرژی روی راکتور توان زیرحوزه را حفظ می‌کند.

**۳. C++ و OpenFOAM:** کرنل سبک `v11/cpp/kinetic_stoichiometry.cpp` کامپایل و تست شد. یک **کاندیدای جدید (کامپایل‌نشده)** از حلگر OpenFOAM در `v11/openfoam/solver` آماده شده که برخلاف v7 شیمی LiCoO₂ دارد و از `pimpleControl` استفاده می‌کند.

**۴. بررسی کیفیت:** ۱۲۴ آزمون Python از v4 تا v11 پاس شده‌اند. علاوه بر آن، مطالعه زمان و شبکه و مقایسه روشن/خاموش شدن منبع EM ثبت شده‌اند.

## راه‌اندازی روی Linux Mint

```bash
sudo apt install g++ make python3 python3-venv
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt pytest
make v11all
```

آزمون جداگانه OpenFOAM فقط بعد از نصب و فعال‌سازی **OpenFOAM.com v2606**:

```bash
make v11openfoam
```

`v11/results/openfoam_runtime/BLOCKED.txt` نشان می‌دهد این بخش در محیط فعلی اجرا نشده است. برای اجرای GitHub Actions، مخزن را با همه فایل‌ها در GitHub آپلود کنید و Workflow مربوط به v11 را اجرا کنید.

## شواهد مهم

- `v11/results/v11_summary.json` — خلاصه قابل پردازش
- `v11/results/integrated_comparison.png` — با و بدون گرمایش القایی
- `v11/results/grid_study.png` — حساسیت به تعداد سلول
- `v11/results/time_step_study.png` — حساسیت به گام زمانی
- `v11/results/all_tests.log` — ۱۲۴ آزمون موفق
- `v11/docs/VERIFICATION_AND_LIMITS_FA.md` — گزارش علمی، واحدها، محدودیت‌ها

**تذکر:** ضرایب سینتیکی از یک مقاله آمده‌اند و به‌شکل تقریب مؤثر در راکتور استفاده شده‌اند؛ این پروژه هنوز اعتبارسنجی صنعتی و OpenFOAM عملیاتی ندارد. حداکثر انحراف فشار محاسباتی حدود ۶٫۷۶٪ است و باید در نسخه بعدی با حل پیوستگی/فشار اصلاح شود.
