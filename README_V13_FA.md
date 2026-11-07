# ReInduct-CFD Pro v13 — راهنمای فارسی

این نسخه، توسعه v12 است و یک **مدل عددی کاهش‌یافته** برای راکتور گاز–جامد و بازیافت LiCoO₂/C ارائه می‌کند. مدل صنعتی، پلاسما–MHD کامل یا حلگر تأیید‌شده OpenFOAM نیست.

## اصلاحات اصلی

- حل تکرارشونده فشار و آنتالپی محسوس گاز در یک گام زمانی؛ کنترل خطای EOS.
- ویژگی‌های حرارتی CO/CO₂/N₂/H₂ بر مبنای Shomate از NIST، با دقت در پیوستگی ضرایب قطعه‌ای؛ H تک‌اتمی ایده‌آل.
- تکرار همگرایی گرمایش ژولی و دمای جامد در هر نوبت به‌روزرسانی EM.
- مقایسه کد C++17 با Python و آزمون‌های مربوط به انرژی، فشار، موازنه اتم‌ها و همگرایی.
- مطالعه ده سناریوی حساسیت؛ **این سناریوها شبیه‌سازی آزمایش مستقل نیستند**.

## اجرای Linux Mint

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt pytest
make v13all
# مطالعه پارامتری کامل جداگانه:
make v13study
```

برای اجرای تنها هسته v13:

```bash
python -m v13.reactor --duration 5 --output v13/results/my_run
python -m v13.study
python -m pytest -q v13/tests
```

برای بررسی **حلگر تاریخی OpenFOAM v11 (غیرمعادل v13)**، جداگانه پس از نصب OpenFOAM.com v2606: `bash v12/openfoam/Allrun`. اجرای `v13/openfoam/Allrun` عمداً با خطا متوقف می‌شود تا برچسب ساختگی Verified تولید نشود.

خروجی‌ها: `v13/results/final_5s/evidence.json`، `v13/results/study/scenarios.csv`، `v13/results/study/conversion_sensitivity.png` و `v13/results/EVIDENCE_MANIFEST.json`.

برای وضعیت دقیق مسائل رفع‌شده و باز، **`v13/docs/ACCEPTANCE_AND_GAPS_FA.md`** را بخوانید.
