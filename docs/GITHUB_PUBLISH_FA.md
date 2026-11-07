# انتشار ReInduct-CFD در GitHub

## روش استاندارد با Git روی Linux Mint

۱. فایل ZIP مخصوص GitHub را استخراج کنید و وارد پوشه `ReInduct-CFD` شوید.
۲. در GitHub یک مخزن **خالی** با نام پیشنهادی `ReInduct-CFD` بسازید (گزینه‌های افزودن README، .gitignore و License را انتخاب نکنید؛ این فایل‌ها موجودند).
۳. از داخل پوشه پروژه دستورات زیر را اجرا کنید:

```bash
make verify

git init -b main
git add .
git status --short
git commit -m "Initial release: ReInduct-CFD v13.1 research prototype"
git remote add origin https://github.com/YOUR_USERNAME/ReInduct-CFD.git
git push -u origin main
```

به جای `YOUR_USERNAME` نام واقعی حساب GitHub خودتان را وارد کنید. در صورت نیاز، `git config --global user.name` و `git config --global user.email` را خودتان تنظیم کنید. هنگام احراز هویت با HTTPS، از credential manager یا token مجاز استفاده کنید؛ رمز یا token را داخل پروژه قرار ندهید.

۴. در زبانه **Actions** مخزن، workflow با عنوان `Numerical verification (not OpenFOAM)` را بررسی کنید. سبز شدن آن فقط نشان می‌دهد تست‌های عددی و C++ در محیط GitHub اجرا شده‌اند؛ به معنی اجرای OpenFOAM نیست.

۵. در بخش About مخزن، متن پیشنهادی زیر مناسب است:

> Reproducible reduced-order reactive gas–solid battery-recycling model with induction-heating feedback, numerical conservation tests and Python/C++17 verification.

موضوعات (Topics) پیشنهادی: `computational-chemistry`, `reactive-transport`, `battery-recycling`, `multiphysics`, `scientific-computing`, `python`, `cpp`.

پیش از انتشار هر داده یا مقاله، مجوز بازنشر و نام نویسندگان را بررسی کنید. نتایج پروژه یک نمونه‌کار عددی‌اند، نه نتایج آزمایشگاهی یا یک Digital Twin اعتبارسنجی‌شده.
