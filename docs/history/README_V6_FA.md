# ReInduct-CFD Pro v6 — نمونه‌کار پژوهشی TU-BRIDGE

**وضعیت: مدل‌های محاسباتی مستقل v4–v6 در محیط محلی آزموده شده‌اند؛ کد OpenFOAM هنوز کامپایل یا اجرا نشده است.**

هدف این نسخه پوشش مستند و صادقانه‌ی مهارت‌های موردنیاز دکتری **DP.MUL.3 (Battery Recycling and Plasma/Inductive Reactor Modeling)** است. پروژه یک **نمونه‌کار آموزشی/تحقیقاتی** است، نه دوقلوی دیجیتال تأییدشده‌ی بازیافت صنعتی.

## دستاوردهای جدید نسخه ۶

۱. **الکترومغناطیس محورمتقارن (Axisymmetric MQS) با Python/SciPy:** مدل پتانسیل برداری حلقوی `A_theta = r * psi`، تحریک حلقوی، جریان‌های گردابی، چگالی توان ژولی `q=0.5 sigma |E|²` و میدان‌های `B_r, B_z`. دستگاه با حل ساختگی تحلیلی (manufactured solution) آزموده شد و همگرایی مرتبه دوم تقریبی دارد. این مدل از سیم‌پیچ کارتزین دو شاخه نسخه ۵ واقع‌گراتر است، اما هندسه سیم‌پیچ واقعی یا پلاسما نیست.

۲. **بستر استوانه‌ای گاز–جامد با FV محافظه‌کار:** حل دوبعدی در مختصات `(z,r)` برای هیدروژن، بخار آب، تبدیل جامد و دمای جداگانه گاز و جامد؛ نیروی محرک واکنش از ترمودینامیک NiO/H₂؛ گرمایش ژولی به دمای جامد منتقل می‌شود؛ رسانایی الکتریکی به دما بازخورد می‌دهد. سرعت محوری از معادله افت فشار **Ergun** به دست می‌آید، نه از حل ناویر–استوکس. مدل بدون لغزش ذرات و بدون واکنش پلاسمایی است.

۳. **ترازهای جرم و انرژی قابل بازسازی:** در هر گام زمانی شار ورودی/خروجی، مصرف H₂، تولید H₂O، تبدیل جامد، تولید گرمای واکنش و توان القایی جمع زده می‌شود. پایستگی در حالت واکنش روشن/خاموش و گرمایش روشن/خاموش آزموده شده است. خطای نزدیک صفر فقط درستی دفترکل عددی را نشان می‌دهد.

۴. **بنچمارک علمی سخت‌گیرانه‌تر:** آزمون ۴-بخشی leave-one-H₂-concentration-out روی ۱۹ داده منتشرشده درباره نرخ تشکیل لایه نیکل در کاهش NiO؛ سه مدل سینتیکی پیش از ارزیابی به‌صورت صریح معرفی شده‌اند. مدل دو‌رژیمی با مرتبه اول، خطای MAE تجمیعی **0.2663 µm/s** دارد، در برابر **0.5646 µm/s** برای مدل آرنیوس ساده. این نرخ رشد سطحی برابر نرخ واکنش حجمی بستر یا بازده بازیافت باتری نیست.

۵. **مدل مواد واقعی باتری:** تراز عنصری و جرمی واکنش کلّی فرض‌شده `4LiCoO2 + 3C -> 2Li2CO3 +4Co+CO2` و سناریوهای مصرف کربن، تبدیل و جداسازی. تبدیل و بازده جداسازی *ورودی فرضی‌اند*، نه مقدار پیش‌بینی‌شده/اندازه‌گیری‌شده.

۶. **C++17:** کد سینتیک مستقل با پارامترهای واحددار و مقایسه عددی با نسخه Python. CI علاوه بر آزمون‌های Python کد C++ را با اخطارهای سخت‌گیرانه کامپایل می‌کند.

۷. **OpenFOAM:** بررسی ساختار فایل‌ها، تعداد سلول‌ها و واحد `qEM` موفق است. گردش‌کار GitHub Actions، نصب OpenFOAM v2606، کامپایل، مش‌بندی و اجرای حلگر را روی Runner تست می‌کند؛ **این گردش‌کار هنوز روی GitHub اجرا نشده است** و هیچ نتیجه OpenFOAM تأییدشده‌ای وجود ندارد.

## نتایج مهم محاسبات، نه آزمایش

- **آزمون‌های خودکار Python/C++ از نسخه ۴ تا ۶:** ۴۳ از ۴۳ موفق.
- **خطای نسبی حل ساختگی محورمتقارن:** 0.0002924 (بدون واحد) در شبکه 40×64.
- **تغییر توان ژولی انتگرالی بین دو شبکه آخر:** حدود 0.0744٪.
- **کیس نماینده ۲ ثانیه‌ای با EM:** بیشینه دمای جامد **942.11 K** و میانگین تبدیل NiO **1.676٪**.
- **بدون گرمایش القایی در شرایط مشابه:** بیشینه دمای جامد **850.04 K**.
- **این مقادیر پیش‌بینی صنعتی نیستند:** هندسه، خواص و سینتیک برای نمایش عددی انتخاب شده‌اند.

## اجرای محلی در Linux Mint

```bash
sudo apt install g++ make python3 python3-venv python3-pip
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
make v6all
# تست تجمیعی
python -m pytest -q v4/tests/test_v4.py v5/tests/test_v5.py v6/tests/test_v6.py
```

نمودارها و JSON خروجی در `v6/results/` قرار می‌گیرند. اجرای `make v6openfoam` به نصب صحیح OpenFOAM v2606 یا Docker نیاز دارد؛ این فرمان تا پایان موفقِ `wmake + blockMesh + checkMesh + solver` نباید موفقیت علمی OpenFOAM را اعلام کند.

## محدودیت‌های مهم برای رزومه

**مجاز:** «Developed and numerically verified an axisymmetric induction-heated gas–solid reactor surrogate with conservative species/energy balances, a C++17 kinetic module, Ergun flow closure and group-heldout NiO kinetic benchmarking; prepared an OpenFOAM v2606 solver candidate (not yet executed).»

**غیرمجاز:** «Validated a plasma-assisted lithium-ion battery recycling plant in OpenFOAM» یا «Demonstrated experimental metal recovery efficiency». داده واقعی برای راکتور صنعتی، مدل کامل پلاسما/MHD و اجرای OpenFOAM هنوز موجود نیست.

## منابع علمی اصلی

- Hidayat et al. 2009, NiO/H₂ measured Ni-layer kinetics, DOI: https://doi.org/10.1007/s11663-008-9212-0
- Review of Li-ion battery recycling, DOI: https://doi.org/10.1039/D3CS00898C
- OpenFOAM v2606 Quickstart: https://doc.openfoam.com/2606/quickstart/
- GitHub Action setup-openfoam: https://github.com/gerlero/setup-openfoam
