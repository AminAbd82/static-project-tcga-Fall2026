# پروژه آمار - پاییز ۱۴۰۵

مسیر ب: **TCGA-BRCA**، اطلاعات بالینی به همراه بیان ژن با RNA-seq.

نسخهٔ تحویل فاز اول بر پایهٔ GDC Data Release 46.0 شامل **۱٬۰۹۴ بیمار دارای هر دو نوع داده**، **۱٬۰۵۴ بیمار کامل در همهٔ متغیرهای اصلی**، **هفت نمودار** و **گزارش فارسی ۱۰ صفحه‌ای** است. گزارش از نتایج محاسبه‌شده تولید می‌شود.

سؤال اولیه: بیان `MKI67` با سن هنگام تشخیص، مرحلهٔ بیماری و بیان `ESR1`، `PGR` و `ERBB2` چه ارتباطی دارد؟ در فاز اول فقط توزیع‌ها و تفاوت‌های توصیفی بررسی می‌شوند.

## فایل‌ها

- `PROJECT_CARD.md`: تعریف سؤال، متغیرها، قواعد انتخاب و محدودیت‌ها
- `data/README.md`: منبع و روش بازسازی داده
- `config/phase1.json`: انتخاب ژن‌ها، واحدها و قواعد پیش از تحلیل
- `src/fetch_data.py`: دریافت داده‌های عمومی و بررسی MD5 فایل‌ها
- `src/analyze_phase1.py`: اتصال، ممیزی کیفیت، جدول‌ها و نمودارها
- `src/build_report.py`: ساخت گزارش فارسی از نتایج واقعی
- `reports/phase1.pdf`: گزارش نهایی فاز اول
- `reports/phase1/`: جدول‌ها، ممیزی و توضیح نتایج
- `figures/phase1/`: نمودارها با عنوان، محور و مقیاس
- `reports/phase1/LEARNING_NOTES_FA.md`: توضیح مفاهیم برای مطالعه و دفاع
- `reports/phase1/REQUIREMENTS_CHECKLIST_FA.md`: تطبیق خروجی‌ها با الزامات درس

## اجرا

نیازمندی: Python 3.10 یا جدیدتر، و `curl` برای دریافت داده.

```bash
python3 -m pip install -r requirements.txt
python3 src/fetch_data.py
python3 src/analyze_phase1.py
python3 src/build_report.py
```

دریافت داده قابل ادامه‌دادن است: فایل‌های کامل موجود دوباره دانلود نمی‌شوند. برای موارد دارای چند فایل RNA، فایل تومور اولیهٔ منتخب مستقیماً دریافت می‌شود؛ فایل خام کامل فقط برای این موارد و نمونه‌های کنترل دریافت خواهد شد. داده‌های خام و محیط نرم‌افزاری داخل GitHub قرار نمی‌گیرند.

## اصل‌های تحلیل

هر بیمار یک ردیف دارد. سن از روز به سال تبدیل می‌شود. مرحلهٔ بیماری متغیر ترتیبی است. واحد بیان ژن FPKM-UQ و تبدیل نمایشی آن `log2(FPKM-UQ + 1)` است. دادهٔ صفر گمشده محسوب نمی‌شود. حذف، گمشدگی و عدم تطابق مستند می‌شوند. تفاوت توصیفی، معناداری آماری یا علیت را اثبات نمی‌کند.

## منابع

- [پروژهٔ TCGA-BRCA](https://portal.gdc.cancer.gov/projects/TCGA-BRCA)
- [روش بیان ژن GDC](https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/)
- [رابط دریافت بیان ژن](https://docs.gdc.cancer.gov/API/Users_Guide/Data_Analysis/)
- [شناسه‌های TCGA](https://docs.gdc.cancer.gov/Encyclopedia/pages/TCGA_Barcode/)
- [پژوهش دربارهٔ چهار ژن منتخب](https://pubmed.ncbi.nlm.nih.gov/40666735/)
