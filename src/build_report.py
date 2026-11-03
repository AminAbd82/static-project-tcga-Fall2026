"""Build a Persian Phase 1 report from the verified, computed results."""
from __future__ import annotations

import json
import math
import sys
from datetime import datetime, timedelta, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if (ROOT / ".deps").exists():
    sys.path.insert(0, str(ROOT / ".deps"))

import arabic_reshaper
from bidi.algorithm import get_display
from reportlab.lib import colors
from reportlab.lib.enums import TA_RIGHT
from reportlab.lib.pagesizes import A4
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import SimpleDocTemplate, Flowable, Table, TableStyle, Image, PageBreak, Spacer

OUT = ROOT / "reports/phase1"
RESULTS = json.loads((OUT / "results.json").read_text())
WIDTH = A4[0] - 88
NAVY = colors.HexColor("#173b4c")
TEAL = colors.HexColor("#256b83")
PALE = colors.HexColor("#eff5f7")
INK = colors.HexColor("#22313a")
pdfmetrics.registerFont(TTFont("Vazir", str(ROOT / "assets/fonts/Vazirmatn-Regular.ttf")))


def shaped(text):
    return get_display(arabic_reshaper.reshape(str(text)), base_dir="R")


class PersianText(Flowable):
    """Wrap logical text before applying Arabic shaping to each displayed line."""
    def __init__(self, text, size=10, color=INK, leading=None, space_after=7):
        super().__init__()
        self.text, self.size, self.color = str(text), size, color
        self.leading = leading or size * 1.75
        self.spaceAfter = space_after

    def wrap(self, available_width, available_height):
        self.width = available_width
        self.lines = []
        for para in self.text.split("\n"):
            line = ""
            for word in para.split():
                test = (line + " " + word).strip()
                if line and pdfmetrics.stringWidth(shaped(test), "Vazir", self.size) > available_width:
                    self.lines.append(line)
                    line = word
                else:
                    line = test
            self.lines.append(line)
        self.height = max(1, len(self.lines)) * self.leading
        return self.width, self.height

    def draw(self):
        self.canv.setFillColor(self.color)
        self.canv.setFont("Vazir", self.size)
        for index, line in enumerate(self.lines):
            self.canv.drawRightString(self.width, self.height - (index + .8) * self.leading, shaped(line))


class SourceLink(Flowable):
    def __init__(self, label, url):
        super().__init__()
        self.label, self.url = label, url
        self.height = 17

    def wrap(self, w, h):
        self.width = w
        return w, self.height

    def draw(self):
        self.canv.setFont("Vazir", 9)
        self.canv.setFillColor(TEAL)
        self.canv.drawString(0, 5, self.label)
        self.canv.linkURL(self.url, (0, 0, self.width, self.height), relative=1)


story = []
narrative = ["# گزارش فاز اول - TCGA-BRCA\n"]


def text(value, size=10, color=INK):
    story.append(PersianText(value, size=size, color=color))
    narrative.append(str(value) + "\n")


def heading(value):
    story.append(Spacer(1, 3))
    text(value, 16, NAVY)
    story.append(Spacer(1, 5))
    narrative.append("\n")


def subheading(value):
    text(value, 11.5, TEAL)


def table(headers, rows, widths, size=9):
    contents = [[PersianText(v, size=size, color=colors.white, leading=size * 1.55, space_after=0) for v in headers]]
    contents += [[PersianText(v, size=size, leading=size * 1.55, space_after=0) for v in row] for row in rows]
    item = Table(contents, colWidths=widths, repeatRows=1, hAlign="CENTER")
    item.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), NAVY), ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, PALE]),
        ("LINEBELOW", (0, 0), (-1, 0), .5, NAVY),
        ("GRID", (0, 1), (-1, -1), .25, colors.HexColor("#d3dfe4")),
        ("LEFTPADDING", (0, 0), (-1, -1), 7), ("RIGHTPADDING", (0, 0), (-1, -1), 7),
        ("TOPPADDING", (0, 0), (-1, -1), 5), ("BOTTOMPADDING", (0, 0), (-1, -1), 5)]))
    story.extend([item, Spacer(1, 9)])
    narrative.append(" | ".join(headers) + "\n")
    narrative.extend(" | ".join(map(str, row)) + "\n" for row in rows)


def figure(index, width=WIDTH, height_limit=235):
    figure = RESULTS["figures"][index]
    from PIL import Image as PILImage
    with PILImage.open(ROOT / figure["file"]) as im:
        w, h = im.size
    scale = min(width / w, height_limit / h)
    image = Image(str(ROOT / figure["file"]), width=w * scale, height=h * scale)
    image.hAlign = "CENTER"
    story.extend([image, Spacer(1, 4)])
    text(f"شکل {index + 1}. {figure['title']}. {figure['caption']}", 8.5)
    narrative.append(f"![{figure['title']}](../../{figure['file']})\n")


def page():
    story.append(PageBreak())


def f(value, decimals=2):
    if value is None:
        return "-"
    return f"{value:,.{decimals}f}" if isinstance(value, (int, float)) else str(value)


def stat(variable, population="multimodal_cohort"):
    return next(r for r in RESULTS["quantitative"] if r["population"] == population and r["variable"] == variable)


def count(step):
    return next(r["n"] for r in RESULTS["flow"] if r["step"] == step)


def categorical(variable, population="multimodal_cohort"):
    return [r for r in RESULTS["categorical"] if r["variable"] == variable and r["population"] == population]


age = stat("age_years")
main_n = RESULTS["n_main"]
received = datetime.fromisoformat(json.loads((OUT / "provenance.json").read_text())["retrieved_at_utc"]).astimezone(timezone(timedelta(hours=3, minutes=30)))
release = RESULTS["gdc_status"]["data_release"]

# Page 1: concise findings and declared scope.
text("پروژهٔ درس آمار | پاییز ۱۴۰۵", 12, TEAL)
heading("فاز اول: شناخت داده و آمار توصیفی")
text("مسیر ب: TCGA-BRCA | اطلاعات بالینی و بیان ژن RNA-seq", 13, NAVY)
text("موضوع: توصیف ویژگی‌های بالینی و بیان ژن‌های منتخب در بیماران سرطان پستان")
text("سؤال اصلی پروژه: بیان MKI67 با سن هنگام تشخیص، مرحلهٔ بیماری و بیان ESR1، PGR و ERBB2 چه ارتباطی دارد؟ این گزارش، بخش توصیفیِ لازم برای بررسی این سؤال را فراهم می‌کند.")
table(["مقدار", "شاخص"], [
    [str(count("project_case_records")), "بیماران ثبت‌شده در پروژه"],
    [str(main_n), "بیماران دارای هر دو نوع داده، پس از انتخاب یک نمونه برای هر بیمار"],
    [str(RESULTS["n_complete"]), "بیماران دارای همهٔ متغیرهای اصلی بدون گمشدگی"],
    [f(age["median"]), "میانهٔ سن هنگام تشخیص، سال"],
    [str(len(RESULTS["figures"])), "تعداد شکل‌های توصیفی"],
    [str(RESULTS["validation_checks"]), "تطبیق عددی موفق بیان ژن با فایل اصلی در موارد کنترل"]], [105, WIDTH - 105])
subheading("دامنهٔ این تحویل")
text("فقط فاز اول اجرا شده است: منشأ داده، واحد مشاهده، قابلیت تعمیم، طبقه‌بندی متغیرها، کنترل کیفیت، جدول‌های فراوانی، مرکز و پراکندگی، مصورسازی و مقایسهٔ توصیفی. هیچ آزمون فرض، مقدار p، برازش رگرسیون یا ادعای علیت در این گزارش وجود ندارد.")
text("منبع: داده‌های عمومی GDC؛ فقط نمونه‌های تومور اولیه در جدول مولکولی اصلی وارد شده‌اند. مقدار صفرِ بیان ژن معتبر است و نقاط دورافتاده صرفاً به‌دلیل بزرگی حذف نشده‌اند.")
text(f"دریافت به وقت تهران در روز {received.day} از ماه {received.month} سال {received.year}، ساعت {received:%H:%M}. نسخهٔ منبع: {release}.", 9)
page()

# Page 2: study design and variable dictionary.
heading("۱. منشأ داده، طراحی مطالعه و متغیرها")
text("TCGA برنامهٔ مشترک NCI و NHGRI است که از سال ۲۰۰۶ آغاز شد و طی ۱۲ سال داده‌های مولکولی و بالینی انواع سرطان را گردآوری کرد [۱]. داده‌های این پروژه از بخش عمومی TCGA-BRCA در GDC دریافت شده‌اند. نسخهٔ انتشار و تاریخ دریافت با زمان تشخیص یا نمونه‌گیری یکسان نیستند.")
text(f"در تشخیص‌های اصلیِ دارای سال معتبر، سال تشخیص از {RESULTS['diagnosis_year_min']} تا {RESULTS['diagnosis_year_max']} است؛ برای {RESULTS['diagnosis_year_available_n']} بیمار این مقدار موجود است. این بازه، بازهٔ سال تشخیص است و نباید به‌جای دورهٔ جمع‌آوری همهٔ نمونه‌ها گزارش شود. تاریخ تقویمی دقیق جمع‌آوری همهٔ نمونه‌ها در استخراج مورد استفاده موجود نیست.")
text("این کار تحلیل ثانویهٔ مشاهده‌ای از پرونده‌های بالینی و اندازه‌گیری‌های مولکولی است. افراد به‌صورت تصادفی برای مداخله تخصیص نیافته‌اند و نمونه، پیمایش تصادفی از جمعیت عمومی نیست. در دسترس‌بودن بافت و دادهٔ باکیفیت می‌تواند ترکیب نمونه را تغییر دهد؛ بنابراین تعمیم به همهٔ بیماران سرطان پستان، به‌خصوص جمعیت ایران، محدود است.")
table(["واحد و نقش", "نوع", "متغیر"], [
    ["کلید اتصال؛ بدون خلاصهٔ عددی", "اسمی، شناسه", "case UUID / بارکد بیمار"],
    ["سال؛ توضیحی", "کمی پیوسته", "سن هنگام تشخیص"],
    ["I تا IV؛ توضیحی و گروه‌بندی", "ترتیبی", "مرحلهٔ پاتولوژیک"],
    ["FPKM-UQ؛ پاسخ", "کمی پیوسته", "بیان MKI67"],
    ["FPKM-UQ؛ توضیحی مولکولی", "کمی پیوسته", "بیان ESR1، PGR، ERBB2"],
    ["برچسب منبع؛ توصیف نمونه", "اسمی", "جنس ثبت‌شده"]], [WIDTH * .40, WIDTH * .23, WIDTH * .37], 8.5)
text("بیان ژن با روش STAR - Counts و نرمال‌سازی FPKM-UQ استفاده شده است [۲،۳]. برای نمایش، تبدیل log2(FPKM-UQ + 1) محاسبه شده؛ تبدیل روی مقدار نرمال‌شده است، نه شمارش خام. بیان RNA با نتیجهٔ آزمایش پروتئینی گیرنده‌ها یکسان فرض نمی‌شود.")
text("انتخاب محدود چهار ژن پیش از محاسبهٔ نتایج ثبت شد. این ژن‌ها در پژوهش‌های بیان نشانگرهای سرطان پستان بررسی شده‌اند [۵]؛ مطالعهٔ حاضر هیچ آستانهٔ تشخیصی یا ادعای بالینی از آن‌ها استخراج نمی‌کند.")
page()

# Page 3: flow and sample selection.
heading("۲. ساختار داده و اتصال در سطح بیمار")
text("ساختار داده سلسله‌مراتبی است: بیمار ← نمونهٔ بافتی ← aliquot ← فایل اندازه‌گیری. یک بیمار می‌تواند چند تشخیص بالینی یا چند فایل RNA داشته باشد. کلید اتصال case UUID است و تطابق بارکد نمونه با بارکد بیمار نیز کنترل شده است [۴].")
labels = {
    "project_case_records": "رکورد بیمار در پروژه", "cases_with_clinical_records": "بیماران دارای رکورد بالینی/جمعیت‌شناختی",
    "clinical_diagnosis_records_before_selection": "رکوردهای تشخیص، پیش از انتخاب تشخیص اصلی",
    "cases_with_unique_primary_diagnosis": "بیماران دارای تشخیص اصلی یکتا",
    "open_star_counts_files_all_tissues": "فایل‌های RNA عمومی، همهٔ انواع بافت",
    "primary_tumor_files_before_patient_selection": "فایل‌های تومور اولیه، پیش از انتخاب یک فایل",
    "non_primary_files_excluded": "فایل‌های بافت طبیعی/متاستاتیک کنارگذاشته‌شده",
    "additional_primary_files_excluded": "فایل‌های تومور اولیهٔ اضافی کنارگذاشته‌شده",
    "unique_primary_rna_patients": "بیماران یکتای مولکولی پس از انتخاب",
    "main_multimodal_patient_cohort": "بیماران جدول اصلی دارای هر دو وجه",
    "complete_all_main_variables": "بیماران کامل در همهٔ متغیرهای اصلی"}
table(["تعداد", "مرحله"], [[str(r["n"]), labels[r["step"]]] for r in RESULTS["flow"] if r["step"] in labels], [65, WIDTH - 65], 8)
text("نمونهٔ انتخابی باید Primary Tumor باشد؛ در صورت چند نمونه، کمترین بارکد نمونه، سپس بارکد aliquot و سپس UUID فایل انتخاب شد. این قاعده ثابت و مستقل از نتایج است و برتری کیفی نمونه را ادعا نمی‌کند. همهٔ کنارگذاری‌ها در sample_file_audit.csv و انتخاب‌ها در selected_files.csv ثبت شده‌اند.", 9)
figure(5, height_limit=165)
page()

# Page 4: clinical descriptive plots.
heading("۳. توصیف و ممیزی بخش بالینی")
text(f"از {count('cases_with_clinical_records')} بیمار دارای رکورد بالینی، {count('cases_with_unique_primary_diagnosis')} بیمار تشخیص اصلی یکتا دارند. فقط همان تشخیص برای سن و مرحله استفاده شد؛ تشخیص قبلی جایگزین تشخیص اصلی نشد. وجود رکورد بالینی به معنی کامل‌بودن همهٔ فیلدها نیست.")
text(f"در نمونهٔ دارای هر دو وجه، سن برای {age['n']} نفر موجود است: میانگین {f(age['mean'])}، میانه {f(age['median'])} و انحراف معیار {f(age['sd'])} سال؛ چارک اول {f(age['q1'])}، چارک سوم {f(age['q3'])} و IQR برابر {f(age['iqr'])} سال است. دامنه از {f(age['min'])} تا {f(age['max'])} سال است.")
figure(0, height_limit=205)
figure(1, height_limit=205)
text("برای سن، میانگین و انحراف معیار خلاصهٔ عمومی‌اند و میانه و IQR نسبت به مقادیر انتهایی مقاوم‌ترند. برای مرحلهٔ ترتیبی، فراوانی و درصد گزارش می‌شود؛ میانگین شمارهٔ مراحل محاسبه نشده است. مخرج درصد مرحله، کل نمونهٔ دارای هر دو وجه است و گمشده‌ها جدا ثبت شده‌اند.", 9)
page()

# Page 5: molecular summaries and scale.
heading("۴. ممیزی مستقل بخش مولکولی")
text("پس از انتخاب یک نمونهٔ تومور اولیه برای هر بیمار، بیان چهار ژن استخراج شد. برای بیمار دارای تنها یک فایل RNA، خروجی سطح بیمارِ رابط GDC بدون ابهام نمونه استفاده شد. برای بیمار دارای چند فایل، مقادیر مستقیماً از فایل تومور منتخب خوانده شد؛ خروجی سطح بیمار به نمونهٔ دلخواه نسبت داده نشد.")
text(f"شش فایل تک‌نمونه‌ای به‌صورت ثابت برای کنترل مستقیم انتخاب شدند؛ {RESULTS['validation_checks']} مقدار ژن با خروجی رابط تطبیق یافت. فایل‌های مستقیم با MD5 رسمی کنترل شدند و تطابق شناسهٔ ژن با نام ژن بررسی شد. جدول و جزئیات کنترل در expression_validation.csv آمده است.")
table(["بیشینه", "کمینه", "IQR", "میانه", "گمشده", "موجود", "ژن"], [
    [f(stat(g + "_fpkm_uq")["max"]), f(stat(g + "_fpkm_uq")["min"]), f(stat(g + "_fpkm_uq")["iqr"]),
     f(stat(g + "_fpkm_uq")["median"]), str(stat(g + "_fpkm_uq")["missing_n"]),
     str(stat(g + "_fpkm_uq")["n"]), g] for g in ["MKI67", "ESR1", "PGR", "ERBB2"]],
    [76, 65, 75, 75, 65, 65, WIDTH - 421], 8)
text("واحد جدول بالا FPKM-UQ است. تعداد موجود و گمشده به نمونهٔ اصلی دارای هر دو نوع داده مربوط‌اند؛ خلاصه‌های مستقل قبل از اتصال نیز در descriptive_statistics.csv ثبت شده‌اند.", 8.5)
figure(2, height_limit=235)
text("تفاوت شکل توزیع روی دو مقیاس، دلیل گزارش هم‌زمان مقیاس اصلی و لگاریتمی است. برای توزیع دارای دنبالهٔ بلند، میانه و IQR معمولاً تفسیر مقاوم‌تری دارند؛ میانگین، انحراف معیار و چارک‌ها نیز برای شفافیت در جدول‌های همراه گزارش شده‌اند.")
page()

# Page 6: all molecular variables and outliers.
heading("۵. توزیع ژن‌ها و نقاط دورافتاده")
figure(4, height_limit=320)
table(["خارج از حصار IQR", "IQR", "انحراف معیار", "میانه", "میانگین", "ژن"], [
    [str(stat(g + "_log2")["outliers_iqr"]), f(stat(g + "_log2")["iqr"]),
     f(stat(g + "_log2")["sd"]), f(stat(g + "_log2")["median"]), f(stat(g + "_log2")["mean"]), g]
    for g in ["MKI67", "ESR1", "PGR", "ERBB2"]], [100, 73, 94, 73, 79, WIDTH - 419], 8)
text("واحد همهٔ خلاصه‌های این جدول log2(FPKM-UQ + 1) است. چارک‌ها با روش درون‌یابی خطی و انحراف معیار نمونه‌ای با مخرج n−1 محاسبه شده‌اند. دورافتاده به معنی خارج‌بودن از بازهٔ Q1−1.5×IQR تا Q3+1.5×IQR است و معادل خطای اندازه‌گیری نیست.")
text("هیچ مشاهده‌ای صرفاً به‌دلیل قرارگرفتن بیرون حصار IQR حذف نشده است. این پرچم فقط برای بررسی توزیع استفاده می‌شود. توزیع ژن‌ها، مقیاس‌های متفاوت و نسبت مقادیر انتهایی نشان می‌دهند که مقایسهٔ مستقیم ارتفاع بیان بین ژن‌های مختلف نیاز به احتیاط دارد؛ سؤال فعلی بر توزیع هر ژن و ارتباط آن در میان بیماران متمرکز است.")
page()

# Page 7: descriptive between-group comparison.
heading("۶. مقایسهٔ توصیفی MKI67 بین مراحل بیماری")
table(["IQR", "میانه", "میانگین", "بیان موجود", "اندازهٔ گروه", "مرحله"], [
    [f(r["iqr"]), f(r["median"]), f(r["mean"]), str(r["n"]), str(r["group_n"]), r["stage"]]
    for r in RESULTS["mki67_by_stage"]], [80, 85, 85, 90, 100, WIDTH - 440], 8.5)
text("مقیاس جدول log2(FPKM-UQ + 1) است. بیمارِ فاقد مرحلهٔ معتبر فقط از این مقایسه کنار گذاشته می‌شود؛ در خلاصهٔ بیان ژنِ کل نمونه باقی می‌ماند. این تفکیک مانع حذف سراسری و پنهان افراد دارای یک فیلد گمشده می‌شود.")
figure(3, height_limit=285)
groups = [r for r in RESULTS["mki67_by_stage"] if r["n"]]
highest, lowest = max(groups, key=lambda r: r["median"]), min(groups, key=lambda r: r["median"])
text(f"در این نمونه، بیشترین میانه در مرحلهٔ {highest['stage']} با مقدار {f(highest['median'])} و کمترین میانه در مرحلهٔ {lowest['stage']} با مقدار {f(lowest['median'])} دیده می‌شود. این مقایسه صرفاً توصیف نمونه است؛ هم‌پوشانی توزیع‌ها و اندازهٔ نامتوازن گروه‌ها باید در تفسیر در نظر گرفته شود.")
smallest = min(groups, key=lambda r: r["n"])
text(f"کوچک‌ترین گروه، مرحلهٔ {smallest['stage']} با {smallest['n']} بیمار دارای بیان موجود است. پایداری خلاصهٔ این گروه از گروه‌های بزرگ‌تر کمتر قابل اتکاست. هیچ آزمون معناداری یا برآورد اثر تعدیل‌شده در این فاز اجرا نشده است.")
page()

# Page 8: missing modalities, variable-level missingness, selection.
heading("۷. گمشدگی و انتخاب نمونهٔ دارای هر دو وجه")
text(f"از کل بیماران پروژه، {count('clinical_only')} نفر فقط وجه بالینی، {count('molecular_only')} نفر فقط وجه مولکولی و {count('neither_modality')} نفر هیچ‌یک از دو وجه واجد شرایط را ندارند؛ {main_n} نفر دارای هر دو وجه‌اند. نبود یک وجه با گمشدگی یک متغیر داخل وجه موجود متفاوت است.")
table(["میانهٔ سن", "سن موجود", "تعداد", "گروه"], [
    [f(r["age_median"]), str(r["age_n"]), str(r["n"]), {"both": "هر دو وجه", "clinical_only": "فقط بالینی", "molecular_only": "فقط مولکولی", "neither": "هیچ‌کدام"}[r["group"]]]
    for r in RESULTS["availability"]], [100, 90, 80, WIDTH - 270], 8.5)
text("توزیع مرحله و جنس نیز به تفکیک «هر دو وجه» و «فقط بالینی» در clinical_selection_comparison.csv گزارش شده است. چون گروه فاقد RNA بسیار کوچک است، حتی تفاوت ظاهری درصدها شواهد قوی برای نبود یا وجود سوگیری انتخاب نیست؛ آزمون آماری روی آن انجام نشده است.")
comparison = RESULTS["selection_comparison"]
both_stage_ii = next(r for r in comparison if r["group"] == "both" and r["variable"] == "stage" and r["category"] == "II")
only_stage_ii = next(r for r in comparison if r["group"] == "clinical_only" and r["variable"] == "stage" and r["category"] == "II")
both_female = next(r for r in comparison if r["group"] == "both" and r["variable"] == "sex_at_birth" and r["category"] == "female")
only_female = next(r for r in comparison if r["group"] == "clinical_only" and r["variable"] == "sex_at_birth" and r["category"] == "female")
text(f"برای نمونه، مرحلهٔ II در گروه دارای هر دو وجه {both_stage_ii['n']} نفر ({f(both_stage_ii['percent'])} درصد) و در گروه فقط بالینی {only_stage_ii['n']} نفر از {only_stage_ii['denominator']} نفر است. نسبت جنس ثبت‌شدهٔ female به‌ترتیب {f(both_female['percent'])} و {f(only_female['percent'])} درصد است. این‌ها تفاوت‌های توصیفی‌اند.", 9)
figure(6, height_limit=235)
text(f"از {main_n} بیمار دارای هر دو نوع داده، {RESULTS['n_complete']} بیمار همهٔ متغیرهای اصلی را بدون گمشدگی دارند. در این فاز خلاصهٔ هر متغیر با دادهٔ موجود همان متغیر محاسبه می‌شود و مخرج آن صریح است. جایگذاری مقدار گمشده انجام نشده و تعداد نمونهٔ کامل جدا گزارش شده است.")
text(f"در دادهٔ اصلی، {RESULTS['age_obfuscated_n']} بیمار دارای پرچم مبهم‌سازی سن هستند. سن دقیق این افراد از دادهٔ عمومی قابل بازیابی نیست؛ خلاصه‌ها سن ثبت‌شده را توصیف می‌کنند. پرچم مربوط در جدول پردازش‌شده نگه داشته شده است.")
page()

# Page 9: decisions, limitations, final question.
heading("۸. تصمیم‌های کیفیت داده و دامنهٔ ادعا")
table(["تصمیم", "مورد بررسی"], [
    ["نامشخص و گزارش‌نشده گمشده هستند؛ مقدار اصلی نیز حفظ می‌شود.", "مرحلهٔ نامشخص"],
    ["فقط تشخیص اصلی یکتا؛ عدم ترکیب فیلدهای تشخیص‌های متفاوت.", "چند تشخیص برای یک بیمار"],
    ["انتخاب ثابت یک فایل تومور اولیه؛ ثبت همهٔ کنارگذاری‌ها.", "چند نمونه یا فایل"],
    ["اعتبار در محدودهٔ صفر تا ۱۲۰ سال؛ مقدار ناممکن به گمشده تبدیل و ثبت می‌شود.", "سن و واحد آن"],
    ["صفر معتبر است؛ مقدار منفی یا نامتناهی نامعتبر است.", "مقدار بیان ژن"],
    ["پرچم‌گذاری با حصار IQR؛ حذف خودکار انجام نمی‌شود.", "مقادیر دورافتاده"],
    ["گزارش تعداد موجود هر متغیر؛ بدون جایگذاری.", "گمشدگی متغیرها"]], [WIDTH * .70, WIDTH * .30], 8.5)
text(f"تعداد مواردی که به‌عنوان مقدار نامعتبر یا مرحلهٔ ناشناخته نیازمند اقدام در quality_issues.csv ثبت شدند: {RESULTS['quality_issue_count']}. نبود اطلاعات به‌خودی‌خود در این شمارش، مقدار نامعتبر محسوب نشده است؛ گمشدگی در جدول‌های مخصوص گزارش می‌شود.")
subheading("محدودیت‌ها")
text("این نمونه، انتخابی و مشاهده‌ای است. وجود بافت مناسب و پوشش داده‌ها، ناهمگنی مراکز و پردازش‌ها، تفاوت ترکیب سلولی تومور، تفاوت ویرایش‌های مرحله‌بندی و نبود اطلاعات برخی بیماران می‌تواند خلاصه‌ها را تحت تأثیر قرار دهد. قاعدهٔ انتخاب نمونه قابل تکرار است اما انتخاب زیستی یگانه‌ای نیست. در این فاز اصلاح اثر مرکز/دستهٔ آزمایش و تحلیل حساسیت اجرا نشده‌اند.")
text("داشتن هر دو وجه به معنی کامل‌بودن همهٔ متغیرها نیست. داده‌ها رابطهٔ علّی یا اعتبار تشخیصی ژن‌ها را ثابت نمی‌کنند. حتی اگر در فازهای بعد نتیجه‌ای معنادار نشود، درستی روش و گزارش محدودیت‌ها مستقل از آن اهمیت دارد.")
subheading("سؤال تثبیت‌شده پس از شناخت داده")
text("در بیماران TCGA-BRCA دارای نمونهٔ تومور اولیهٔ RNA-seq و اطلاعات بالینی، بیان MKI67 چگونه توزیع شده و از نظر توصیفی میان مراحل پاتولوژیک بیماری چه تفاوتی دارد؟ بررسی ارتباط رسمی آن با سن و سایر ژن‌های از پیش انتخاب‌شده، ادامهٔ سؤال پروژه برای فازهای بعد خواهد بود.")
text("دادهٔ حاضر برای بخش توصیفی سؤال مناسب است، چون اندازه‌گیری کمی بیان ژن، سن و گروه‌بندی مرحله در سطح بیمار قابل اتصال‌اند. ناقص‌بودن مرحله در برخی افراد و کوچک‌بودن بعضی گروه‌ها محدودیت این تناسب است.")
page()

# Page 10: reproducibility and references.
heading("۹. بازتولید، فایل‌های همراه و منابع")
text("ترتیب بازسازی: ابتدا وابستگی‌های requirements.txt نصب می‌شوند؛ سپس src/fetch_data.py داده‌های عمومی را دریافت می‌کند، src/analyze_phase1.py جدول‌ها و شکل‌ها را می‌سازد و src/build_report.py همین گزارش را از نتایج محاسبه‌شده تولید می‌کند. نتایج عددی به‌صورت دستی در گزارش وارد نشده‌اند.")
text("فایل‌های خام در data/raw نگه داشته می‌شوند و در GitHub یا بستهٔ تحویل قرار نمی‌گیرند. پروندهٔ provenance.json نسخه و تاریخ دریافت، source_fingerprints.csv اثر انگشت فایل‌های خام و selected_files.csv شناسه و MD5 منبع را ثبت می‌کنند. یکسان‌بودن نسخهٔ انتشار و اثر انگشت فایل‌ها برای بازتولید دقیق مهم است.")
table(["محتوا", "فایل همراه"], [
    ["جریان نمونه و تعدادها", "data_flow.csv"],
    ["تصمیم برای تمام فایل‌ها و تشخیص‌ها", "sample_file_audit.csv / diagnosis_audit.csv"],
    ["خلاصهٔ کمی قبل و بعد از اتصال", "descriptive_statistics.csv"],
    ["فراوانی طبقه‌ها و مقایسهٔ گروه‌ها", "categorical_frequencies.csv / mki67_by_stage.csv"],
    ["الگوی وجود داده و مقایسهٔ نمونه‌ها", "modality_availability.csv / clinical_selection_comparison.csv"],
    ["کنترل عددی و مقادیر نامعتبر", "expression_validation.csv / quality_issues.csv"]], [WIDTH * .45, WIDTH * .55], 8)
subheading("منابع")
references = [
    ("[1] NCI - The Cancer Genome Atlas", "https://www.cancer.gov/ccg/research/genome-sequencing/tcga"),
    ("[2] GDC - mRNA Analysis Pipeline", "https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/"),
    ("[3] GDC API - Gene Expression Values", "https://docs.gdc.cancer.gov/API/Users_Guide/Data_Analysis/"),
    ("[4] GDC - TCGA Barcode", "https://docs.gdc.cancer.gov/Encyclopedia/pages/TCGA_Barcode/"),
    ("[5] ESR1, PGR, ERBB2 and MKI67 mRNA study (2025)", "https://pubmed.ncbi.nlm.nih.gov/40666735/"),
    ("TCGA-BRCA project and public data", "https://portal.gdc.cancer.gov/projects/TCGA-BRCA")]
for label, url in references:
    story.append(SourceLink(label, url))
    narrative.append(f"- [{label}]({url})\n")
text("قدردانی: این نتایج بر داده‌های تولیدشده توسط شبکهٔ پژوهشی TCGA متکی‌اند. از اهداکنندگان نمونه و گروه‌های گردآورندهٔ داده قدردانی می‌شود.", 9)
text("شفافیت آماده‌سازی: کد، محاسبات و پیش‌نویس گزارش با کمک دستیار هوش مصنوعی آماده شده‌اند؛ بازبینی علمی و توانایی توضیح و دفاع از تصمیم‌ها بر عهدهٔ دانشجو است.", 8.5)


def decorate(canvas, doc):
    canvas.saveState()
    w, h = A4
    canvas.setStrokeColor(colors.HexColor("#d4e1e7"))
    canvas.line(44, 40, w - 44, 40)
    canvas.setFillColor(TEAL)
    canvas.setFont("Vazir", 8)
    canvas.drawRightString(w - 44, 27, shaped("فاز اول | آمار توصیفی | TCGA-BRCA"))
    canvas.drawString(44, 27, str(doc.page))
    canvas.restoreState()


def main():
    path = ROOT / "reports/phase1.pdf"
    doc = SimpleDocTemplate(str(path), pagesize=A4, leftMargin=44, rightMargin=44,
                            topMargin=43, bottomMargin=51, title="TCGA-BRCA - Phase 1 Descriptive Statistics",
                            author="Statistics course project", pageCompression=1)
    doc.build(story, onFirstPage=decorate, onLaterPages=decorate)
    (OUT / "report_fa.md").write_text("\n".join(narrative), encoding="utf-8")
    print(f"Created: {path}")


if __name__ == "__main__":
    main()
