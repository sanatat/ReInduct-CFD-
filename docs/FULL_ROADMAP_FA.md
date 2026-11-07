# ReInduct-CFD Pro — نقشه پوشش هفت شکاف TU-BRIDGE

## وضعیت شفاف (۵ نوامبر ۲۰۲۶)

| شکاف | شواهد قابل اجرا در این بسته | ادامه لازم |
|---|---|---|
| C++ | `src/reactor1d.cpp`، آزمایش‌های موازنه جرم | توسعه سورس‌های CFD در C++ و مستندسازی API |
| Reactive CFD | شبیه‌ساز انتقال واکنشی یک‌بعدی C++؛ **CFD کامل نیست** | حل واقعی Navier–Stokes + گونه‌ها + انرژی در OpenFOAM |
| Multiphase/multiregion | مدل گاز/جامد با دو دمای مستقل (reduced-order) | regionهای مستقل، boundary mapping، سطح مشترک |
| High temperature/induction | حرارت القایی حل‌شده با EM و انتقال گرما در بستر صفحه‌ای | خواص وابسته به T، اتلاف تابشی، هندسه واقعی کویل |
| Electromagnetics/MHD | حل عددی harmonic magnetic diffusion و جریان گردابی، تلفات ژول، مقایسه با جواب تحلیلی | Lorentz force + flow feedback، میدان پلاسما، اعتبارسنجی |
| OpenFOAM | case اولیه دوبعدی incompressible با فایل‌های mesh/BC/solver | روی سیستم دارای OpenFOAM اجرا/اصلاح و به reactive CFD ارتقا یابد |
| Thermodynamics/kinetics | Arrhenius illustrative + bookkeeping reaction | داده واقعی ΔG°(T)، K(T)، نرخ‌های literature و uncertainty |

## تعریف نهایی پژوهش

**Coupled simulation of a gas–solid reactive battery recycling surrogate under induction heating:**
1. گاز H2 و H2O در بستر نماینده NiO/Ni با انتقال گونه و واکنش سطحی.
2. تخلخل/افت فشار (Ergun/Darcy–Forchheimer) و جریان وابسته به دما.
3. حل گذرای انرژی جامد و سیال با تبادل حرارتی دوطرفه.
4. حل مستقل میدان هارمونیک با هدایت الکتریکی در فاز جامد، و نگاشت q=|J|²/(2σ) به مش حرارتی.
5. فاز بعدی: حل مغناطیس‌هیدرودینامیک با J×B و واکنش‌های رادیکالی مبتنی بر مقالات معتبر.

این پروژه یک **surrogate** است؛ باتری واقعی مخلوطی از NMC/LFP، کربن، binder و ناخالصی است و NiO تنها جایگزین اثبات روش است.

## مراحل اجرایی

### M1 — CFD پایه و مش‌بندی
- اجرای `openfoam/cavity2D` با نسخه‌ای که `icoFoam` دارد.
- ثبت OpenFOAM version، residualها، `checkMesh`، Re، Courant و field visualization.
- تکرار 30x30، 60x60، 120x120 و مقایسه U_center. کد اجرای خودکار همه شبکه‌ها هنوز ساخته نشده است.

### M2 — Reacting CFD واقعی
- از tutorial رسمی reacting solver همان distribution کپی کن؛ انتخاب solver به نسخه وابسته است.
- واکنش مدل را به source terms گونه‌های H2,H2O و تغییر موجودی NiO اضافه کن.
- اگر solid جداست، صحت معادله گونه و شرط عدم نفوذ روی سطح مشترک را چک کن.
- تست‌های خاموش بودن واکنش، حفظ اتم‌های H و O، حساسیت نرخ به T و مسیریابی گونه‌ها.

### M3 — Multiregion
- گاز و solid porous region، کپل انرژی با شار مرزی برابر و مخالف.
- شاخص اعتبارسنجی: |Qgas+Qsolid| نسبت به Qinterface کم‌تر از ۱٪ پس از همگرایی عددی.
- تست conduction-only با جواب تحلیلی steady-state.

### M4 — Induction coupling
- کد `em/induction1d.py` یک مسئله صفحه‌ای خطی را حل و نسبت به جواب `cosh` می‌سنجد.
- `em/coupled_thermal.py` با خواص ثابت، توزیع q را به معادله انرژی وارد می‌کند.
- مرحله تکمیلی: استفاده از Elmer FEM برای هندسه کویل سه‌بعدی و واردکردن q به OpenFOAM.

### M5 — Thermodynamics and kinetics
- جمع‌آوری ΔH°, ΔS°, Cp(T), فعالیت‌ها، واکنش‌های برگشت‌پذیر از مقالات و پایگاه‌های معتبر با ارجاع.
- برای NiO+H2 ↔ Ni+H2O، داده‌های تعادلی قابل‌ردیابی، فیت سینتیک غیرایزوترمال و بازه عدم‌قطعیت ضروری است.
- لازم است قبل از استفاده در مقاله/رزومه، نرخ‌های illustrative کد اولیه جایگزین شوند.

### M6 — Plasma and MHD (advanced)
- در محیط plasma پارامترهای رسانایی و یونش از دما/ترکیب تعیین می‌شوند؛ مدل harmonic slab برای plasma کافی نیست.
- پیاده‌سازی سیستم الکترومغناطیس + گاز، نیروی J×B و گرمایش ژول با موازنه انرژی.
- در صورت عدم دسترسی به داده واقعی، فقط benchmarkهای معتبر توصیف شوند نه «پلاسما اعتبارسنجی‌شده».

## معیارهای کیفیت برای رزومه
- ویدئو یا تصویر solver output + شماره commit + نسخه برنامه.
- مقایسه mesh و time-step و residual history؛ test report با PASS/FAIL.
- منابع و citation DOI مربوط به kinetics و material properties.
- نمودارها با محور و واحد، گزارش sensitivity و limitation.
- کدهای تحت گیت و دستور اجرا روی Linux با حداقل منابع.

## متن قابل دفاع CV (پس از بررسی شخصی اجرای کدها)
"Developed and tested a C++ finite-volume surrogate of transient gas–solid redox conversion, with conservation checks and mesh sensitivity analysis; implemented a harmonic magnetic-diffusion benchmark with Joule-heating-to-thermal coupling and analytic verification (Python/NumPy). Prepared a 2-D OpenFOAM incompressible-flow case for subsequent reactive-flow extension."

**ادعا نکن:** "validated industrial plasma reactor", "completed OpenFOAM reactive CFD", "full MHD coupling", or "experimental battery recycling".
