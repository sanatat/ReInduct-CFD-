# ReInduct-CFD Pro — Enhanced portfolio release

این بسته نسخه توسعه‌یافته پروژه اولیه است. دو جزء **جدید و اجراشده** عبارت‌اند از حل عددی میدان مغناطیسی در ورق رسانا و انتقال یک‌طرفه گرمایش القایی به معادله انرژی. کد C++ مدل واکنشی نسخه قبل و تست‌های آن حفظ شده‌اند. یک کیس OpenFOAM پایه نیز آماده شده اما در محیط موجود به علت نصب نبودن OpenFOAM اجرا نشده است.

## اجرای نسخه سبک روی Linux Mint

```bash
sudo apt install g++ make python3 python3-numpy python3-matplotlib
make all test emtest
make run
make em
make coupled
make plots
```

`make openfoam` فقط وقتی OpenFOAM سازگار نصب و environment آن source شده باشد قابل اجرا است.

## شواهد تکمیل
- `src/reactor1d.cpp` — C++17 redox surrogate + mass conservation.
- `em/induction1d.py` — Maxwell-derived harmonic slab benchmark + exact-solution verification.
- `em/coupled_thermal.py` — numerical Joule heat into thermal diffusion.
- `tests/test_induction.py`, `tests/test_coupling.py` — refinement, analytic comparison, zero-source and quadratic H0 tests.
- `openfoam/cavity2D` — non-reactive 2D flow starter (untested here).
- `docs/FULL_ROADMAP_FA.md` — step-by-step gap coverage for TU-BRIDGE.

**یادآوری مهم:** کد هنوز CFD واکنشی کامل، مدل چندفازی OpenFOAM، MHD جریان-میدان و شیمی پلاسما را حل نکرده است. این‌ها هدف مراحل بعدی‌اند و نباید در CV به عنوان تجربه تکمیل‌شده مطرح شوند.
