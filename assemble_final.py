from pathlib import Path
import sys, io, json, re, hashlib, zipfile, shutil
BASE=Path(__file__).resolve().parent
ROOT=BASE
SRC=BASE
import pymupdf as fitz
from reportlab.pdfgen import canvas
from reportlab.lib.utils import simpleSplit

OUT=BASE/'build'
OUT.mkdir(exist_ok=True,parents=True)
doc=fitz.open(SRC/'manuscript.pdf')
title='Logical Slab Access and Local Processing with Three-Dimensional Layer Codes at Optimal Distance'
rows=[]; seen=set()
for level,t,p in doc.get_toc():
    if level!=2 or t in ('Abstract','Acknowledgments','Acknowledgment of AI assistance') or t in seen: continue
    seen.add(t);rows.append((t,p))
buf=io.BytesIO(); c=canvas.Canvas(buf,pagesize=(612,792))
c.setTitle(title);c.setAuthor('Aryaman Katoch')
c.setFont('Times-Roman',10);c.drawCentredString(306,733,'IRIS RESEARCH PAPER')
y=684;c.setFont('Times-Bold',17)
for line in simpleSplit(title,'Times-Bold',17,468):
    c.drawCentredString(306,y,line);y-=23
c.setFont('Times-Roman',12);c.drawCentredString(306,y-15,'Aryaman Katoch')
y-=65;c.setFont('Times-Bold',13);c.drawString(72,y,'Contents')
c.setLineWidth(.4);c.line(72,y-11,540,y-11);y-=33
links=[]
for t,p in rows:
    c.setFont('Times-Roman',10.5);c.drawString(72,y,t);c.drawRightString(540,y,str(p))
    links.append((fitz.Rect(70,792-y-11,542,792-y+3),p))
    y-=21
assert y>125,(y,len(rows))
c.setFont('Times-Italic',9)
c.drawString(72,112,'Page numbers refer to the numbered manuscript body that follows.')
c.drawString(72,98,'Proofs, assumptions, and numerical methods are included in this volume.')
c.setFont('Times-Roman',10);c.drawCentredString(306,45,'i');c.showPage();c.save()
cover=fitz.open(stream=buf.getvalue(),filetype='pdf')
result=fitz.open();result.insert_pdf(cover,links=False);result.insert_pdf(doc,links=False)
names=doc.resolve_names(); nlinks=0
for i,page in enumerate(doc):
    for old in page.get_links():
        base={'from':old['from']}
        if old.get('nameddest'):
            d=names[old['nameddest']]; p=d['page']
            pt=fitz.Point(*d.get('to',(0,doc[p].rect.height)))*doc[p].transformation_matrix
            new=dict(base,kind=fitz.LINK_GOTO,page=p+1,to=pt,zoom=0)
        elif old.get('uri'):
            new=dict(base,kind=fitz.LINK_URI,uri=old['uri'])
        elif old['kind']==fitz.LINK_GOTO:
            new=dict(base,kind=fitz.LINK_GOTO,page=old['page']+1,to=old.get('to',fitz.Point(0,0)))
        else:raise ValueError(old)
        result[i+1].insert_link(new);nlinks+=1
for rect,p in links:result[0].insert_link(dict(kind=fitz.LINK_GOTO,**{'from':rect},page=p,to=fitz.Point(0,0)))
toc=[[1,'Title and contents',1]];seen=set()
for level,t,p in doc.get_toc():
    if (level,t,p) in seen:continue
    seen.add((level,t,p));toc.append([level,t,p+1])
result.set_toc(toc)
result.set_page_labels([{'startpage':0,'prefix':'','style':'r','firstpagenum':1},
                        {'startpage':1,'prefix':'','style':'D','firstpagenum':1}])
result.set_metadata({'title':title,'author':'Aryaman Katoch','creator':'LaTeX; revised manuscript assembly'})
final=OUT/'IRIS_MANUSCRIPT_FINAL.pdf'
result.save(final,garbage=4,deflate=True);result.close()
with fitz.open(final) as f:
    assert len(f)==len(doc)+1
    text='\n'.join(p.get_text() for p in f)
    assert '??' not in text
    assert 'Conditional universal architecture' in text
    for page in f:
        for l in page.get_links():
            if l['kind']==fitz.LINK_GOTO:assert 0<=l['page']<len(f)
            else:assert l['kind']==fitz.LINK_URI,l
    for i,p in enumerate(doc):
        assert p.get_pixmap().samples==f[i+1].get_pixmap().samples
    (BASE/'final_text.txt').write_text(text,encoding='utf-8')
    report={'pages':len(f),'body_pages':len(doc),'body_links':nlinks,'contents_links':len(links),
            'all_internal_links_valid':True,'body_render_identical':True,
            'sha256':hashlib.sha256(final.read_bytes()).hexdigest()}
    (BASE/'pdf_validation.json').write_text(json.dumps(report,indent=2))
print(json.dumps(report,indent=2))
