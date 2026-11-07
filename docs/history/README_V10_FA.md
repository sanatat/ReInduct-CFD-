# ReInduct-CFD Pro v10 — نسخه تقویت‌شده و قابل ممیزی

**تاریخ:** ۶ نوامبر ۲۰۲۶. **هدف علمی:** نمونه‌کار مرتبط با پروژه TU-BRIDGE DP.MUL.3 در زمینه مدل‌سازی راکتور بازیافت حرارتی باتری. این کد **Digital Twin صنعتیِ اعتبارسنجی‌شده نیست**.

## پیشرفت‌های اجراشده در نسخه ۱۰

1. **انتقال فضایی پلاسما:** یک مدل **یک‌بعدی و بی‌بُعد** شامل H₂، H، H₂⁺، الکترون، دو مخزن انرژی و پتانسیل پواسون. معادلات انتقال رانش–نفوذ با حجم محدود حل می‌شوند. تولید/بازترکیب ذرات باردار بار کل را حفظ می‌کند؛ موازنه هیدروژن و انرژی شامل انرژی ذخیره‌شده در میدان الکتروستاتیک آزموده می‌شود. **این مدل یک تخلیه پلاسمایی کالیبره‌شده، شبیه‌سازی RF یا راکتور واقعی نیست.**
2. **MHD و حرارت واقعی‌تر:** حلگر مستقل جریان Hartmann دارای پارامترهای بعددار است و میدان سرعت عددی را با حل تحلیلی مقایسه می‌کند. مقاومت مغناطیسی، اتلاف ویسکوز، گرمایش ژولی و افزایش دما به شکل منسجم محاسبه می‌شوند.
3. **سینتیک واقعی‌تر LiCoO₂/C:** ده ضریب سرعت گزارش‌شده در جدول‌های II و III مقاله Nuraeni و همکاران (۲۰۲۳) مستقیماً وارد برنامه شده‌اند. انرژی فعال‌سازی دو مرحله، و آزمون Leave-One-Temperature-Out روی *ضرایب گزارش‌شده در همان مقاله* محاسبه می‌شوند. آزمون دوم هنوز اعتبارسنجی خارجی راکتور نیست.
4. **سنجش مطالعه مستقل:** داده‌های بازیابی Li₂CO₃ و CoO از مقاله She و همکاران (۲۰۲۲) جداگانه مستند شده‌اند؛ این مقاله در هوا انجام شده، درحالی‌که مطالعه اصلی سینتیک در آرگون بوده است. اختلاف CoO و Co و تفاوت «درصد بازیابی» با «کسر تبدیل» عمداً شناسایی و ثبت شده‌اند؛ عددسازی نکرده‌ایم.
5. **OpenFOAM:** مسیر امن برای تلاش کامپایل و اجرای **واقعی** v2606 با لاگ‌های wmake، blockMesh، checkMesh، exportCellGeometry، نگاشت qEM و حلگر مستقل آماده شده است. پروژه اصلی دست‌نخورده می‌ماند. در این محیط `wmake` موجود نیست؛ **OpenFOAM هنوز اجرا نشده و موفقیت آن ادعا نمی‌شود**.

## نتایج محاسباتی این اجرا

- **۱۰۸ تست Python موفق** (شامل ۱۹ تست جدید نسخه ۱۰)، به‌علاوه کامپایل و اجرای یک آزمون C++ مستقل؛ این آزمون C++ مربوط به خود OpenFOAM نیست.
- خطای بیشینه موازنه اتم هیدروژن در مدل فضایی: تقریباً `1.33e-15` واحد نرمال‌شده.
- خطای بار خالص یکپارچه: حدود `6.11e-18`؛ خطای موازنه انرژی شامل میدان: حدود `2.41e-12`.
- خطای نسبی L2 سرعت Hartmann در شبکه ۹۷ نقطه‌ای: حدود `6.98e-5`؛ در شبکه ۱۰۱ نقطه‌ای: `6.43e-5`.
- انرژی فعال‌سازی از داده‌های منتشرشده: `121.17 kJ/mol` برای مرحله GB و `95.28 kJ/mol` برای مرحله A3.
- OpenFOAM: تلاش اجرا با خروج کد ۲ و پیام `missing executable wmake`؛ بدون خروجی CFD واقعی.

## اجرا

پس از بازکردن بسته و ورود به پوشه پروژه:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
make v10all
```

برای تلاش اجرای واقعی OpenFOAM (پس از نصب و فعال‌سازی v2606):

```bash
make v10openfoam
```

برای اجرای CI، پوشه `.github/workflows` را همراه با تمام فایل‌های پروژه در یک مخزن GitHub قرار دهید؛ سپس از بخش Actions، گردش‌کار v10 را اجرا و لاگ‌ها را بررسی کنید. **ساخته‌شدن فایل Workflow به معنای اجرای آن نیست.**

## مستندات و شواهد

- `v10/results/evidence_manifest.json` — تعداد تست‌ها، شاخص‌ها، وضعیت OpenFOAM و SHA-256 فایل‌ها.
- `v10/results/plasma1d.json` و `plasma1d.npz` — مدل فضایی پلاسما و میدان.
- `v10/results/hartmann.json` و `hartmann.npz` — جریان و انرژی.
- `v10/results/published_rate_table.csv` — داده‌های دقیق نرخ و منشأ.
- `v10/results/rate_validation.json` — برازش و آزمون دماهای کنارگذاشته‌شده.
- `v10/results/external_phase_audit.json` — سنجش کیفی بین دو پژوهش و اختلاف‌ها.
- `v10/results/openfoam_attempt.log` — تلاش واقعی و دلیل توقف.
- `v10/docs/GAPS_AND_ACCEPTANCE_FA.md` — شکاف‌های هنوز بسته‌نشده و معیارهای دقیق تکمیل.

## منابع اصلی

- Nuraeni et al., *Metallurgical and Materials Transactions B*, 2023. DOI: https://doi.org/10.1007/s11663-022-02712-1
- She et al., *Journal of Chemical Research*, 2022. DOI: https://doi.org/10.1177/17475198211066533
- OpenFOAM v2606, official: https://www.openfoam.com/current-release

## ادعای درست برای رزومه

> Developed and tested coupled numerical verification prototypes for gas–solid recycling models, including conservative 1D charge/species/energy transport, a benchmarked Hartmann MHD–thermal solver, and literature-grounded LiCoO₂/C kinetic analysis with documented limitations. Prepared an OpenFOAM v2606 runtime-validation pipeline; successful OpenFOAM execution remains pending.
