"""Draws docs/workflow.pdf and docs/workflow.png (HuntingTool analysis workflow). Run from the repository root: python docs/make_workflow.py"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from matplotlib.patches import FancyBboxPatch, Rectangle, Polygon, Circle
import sys
OUT = sys.argv[1] if len(sys.argv) > 1 else 'docs/workflow'
FONT = 'Inter' if any('Inter' == f.name for f in fm.fontManager.ttflist) else 'DejaVu Sans'
plt.rcParams.update({'font.family': FONT, 'pdf.fonttype': 42, 'svg.fonttype': 'none'})
INK='#1b2430'; MUT='#566170'; LINE='#9aa5b1'; HEAD='#2c3e55'; PALE='#f3f5f8'
BLUE='#2f6fb3'; ORG='#d98324'; GRN='#3b9468'; GRY='#b5bdc7'; RED='#b23b3b'; PUR='#7a5aa6'
W,H=100,56.25
fig=plt.figure(figsize=(13.33,7.5),dpi=200); ax=fig.add_axes([0,0,1,1]); ax.set_xlim(0,W); ax.set_ylim(0,H); ax.axis('off')
def T(x,y,s,fs=7.6,c=INK,w='normal',ha='left',va='center',**k): return ax.text(x,y,s,fontsize=fs,color=c,fontweight=w,ha=ha,va=va,zorder=6,**k)
def box(x,y,w,h,num,title):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0,rounding_size=0.6',fc='white',ec=LINE,lw=1.0,zorder=1))
    ax.add_patch(FancyBboxPatch((x,y+h-3.6),w,3.6,boxstyle='round,pad=0,rounding_size=0.6',fc=HEAD,ec=HEAD,lw=1.0,zorder=2))
    ax.add_patch(Rectangle((x,y+h-3.6),w,1.6,fc=HEAD,ec=HEAD,lw=0,zorder=2))
    ax.add_patch(Circle((x+2.0,y+h-1.8),1.15,fc='white',ec='none',zorder=3)); T(x+2.0,y+h-1.85,num,8.6,HEAD,'bold','center')
    T(x+4.0,y+h-1.8,title,9.0,'white','bold')
def arrow(x0,y0,x1,y1,c=MUT,lw=1.6,head=True):
    ax.annotate('',xy=(x1,y1),xytext=(x0,y0),arrowprops=dict(arrowstyle='-|>' if head else '-',color=c,lw=lw,shrinkA=0,shrinkB=0,mutation_scale=13),zorder=5)
def hline(x0,x1,y,c=LINE): ax.plot([x0,x1],[y,y],lw=0.7,color=c,zorder=4)
def seg(x,y,w,h,c,lab=None,fs=6.6):
    ax.add_patch(Rectangle((x,y),w,h,fc=c,ec='white',lw=1.0,zorder=4))
    if lab: T(x+w/2,y+h/2,lab,fs,'white','bold','center')
# geometry
bw=21.6; gap=3.5; x0=1.55
X=[x0+i*(bw+gap) for i in range(4)]
Y1,H1=29.6,25.6; Y2,H2=1.5,22.5
# ------------------------------------------------ row 1
# 1 Input and filtering
x=X[0]; box(x,Y1,bw,H1,'1','Input and filtering')
T(x+1.2,Y1+H1-6.0,'FASTQ or FASTQ.GZ files, or a MinKNOW',7.6); T(x+1.2,Y1+H1-7.5,'fastq_pass folder (one sample per barcode)',7.6)
hline(x+1.2,x+bw-1.2,Y1+H1-9.0)
T(x+1.2,Y1+H1-10.6,'Read length ≥ 100 bp (user-defined)',7.6); T(x+1.2,Y1+H1-12.1,'Optional minimum mean read quality (Q)',7.6)
ys=Y1+H1-14.4
reads=[(0.0,15.5,1),(2.0,15.0,1),(0.0,7.5,0),(3.0,15.5,1),(0.0,5.5,0),(1.0,16.0,1)]
for i,(a_,b_,keep) in enumerate(reads):
    yy=ys-i*1.55; ax.add_patch(Rectangle((x+1.6+a_,yy-0.3),b_,0.6,fc=BLUE if keep else GRY,ec='none',zorder=4))
ax.add_patch(Rectangle((x+1.6,Y1+2.0),1.4,0.8,fc=BLUE,ec='none',zorder=4)); T(x+3.4,Y1+2.4,'retained',6.6,MUT)
ax.add_patch(Rectangle((x+9.6,Y1+2.0),1.4,0.8,fc=GRY,ec='none',zorder=4)); T(x+11.4,Y1+2.4,'discarded',6.6,MUT)
# 2 Orientation
x=X[1]; box(x,Y1,bw,H1,'2','Orientation')
T(x+1.2,Y1+H1-6.0,'Each read is oriented on the HTT sense',7.6); T(x+1.2,Y1+H1-7.5,'strand: both strands are analysed',7.6)
yf=Y1+H1-11.4
T(x+1.2,yf+1.9,'Forward-strand reads  (CAG)n',7.0,BLUE,'bold'); arrow(x+1.4,yf,x+bw-1.4,yf,BLUE,2.0)
yr=yf-5.2
T(x+1.2,yr+1.9,'Reverse-strand reads  (CTG)n',7.0,PUR,'bold'); arrow(x+bw-1.4,yr,x+1.4,yr,PUR,2.0)
arrow(x+bw/2,yr-1.6,x+bw/2,yr-3.6,MUT,1.6)
ax.add_patch(FancyBboxPatch((x+1.4,yr-7.0),bw-2.8,2.9,boxstyle='round,pad=0,rounding_size=0.4',fc=PALE,ec=LINE,lw=0.8,zorder=3))
T(x+bw/2,yr-5.55,'reverse-complement → sense orientation',7.0,INK,'normal','center')
T(x+1.2,Y1+1.5,'Reverse-strand reads: ≈38% of parsed reads',6.6,MUT,style='italic')
# 3 Anchoring
x=X[2]; box(x,Y1,bw,H1,'3','Anchoring')
T(x+1.2,Y1+H1-6.0,'Each read is anchored on three sequences',7.6)
bx=x+1.3; by=Y1+H1-12.0
parts=[(2.8,GRY,'flank'),(5.8,BLUE,'(CAG)n'),(3.6,ORG,'adj.'),(4.6,GRN,'(CCG)n'),(2.4,GRY,'CCT')]
cx=bx; pos=[]
for ww,c,l in parts: seg(cx,by,ww,2.2,c,l); pos.append(cx); cx+=ww
pos.append(cx)
def mark(xp,n):
    ax.add_patch(Polygon([[xp-0.55,by+3.0],[xp+0.55,by+3.0],[xp,by+2.45]],fc=RED,ec='none',zorder=5))
    ax.add_patch(Circle((xp,by+4.0),0.9,fc='white',ec=RED,lw=0.9,zorder=5)); T(xp,by+3.95,n,6.6,RED,'bold','center')
mark(pos[1],'1'); mark(pos[3],'2'); mark(pos[4],'3')
for i,(n,l) in enumerate([('1','5′ flank CTCAAGTCCTTC (1 mismatch tolerated)'),('2','start of the CCG repeat'),('3','CCT anchor')]):
    yy=by-2.9-i*2.1; ax.add_patch(Circle((x+1.9,yy),0.95,fc='white',ec=RED,lw=0.9,zorder=5)); T(x+1.9,yy-0.05,n,6.6,RED,'bold','center'); T(x+3.5,yy,l,7.0,INK)
T(x+1.2,Y1+4.1,'Anchors 1–2 delimit the CAG tract and the',6.9,MUT); T(x+1.2,Y1+2.7,'adjacent sequence (adj.); CAG count = tract / 3',6.9,MUT)
# 4 Configuration
x=X[3]; box(x,Y1,bw,H1,'4','Configuration')
T(x+1.2,Y1+H1-5.9,'Adjacent sequence: exact match after ≥2 CAG codons',7.0)
cfg=[('canonical','CAACAGCCGCCA'),('CAA-loss','CCGCCA'),('CCA-loss','CAACAG'),('CAACAGCCGCCA-loss','none'),('CAACAG-dup (G2)','CAACAGCAACAGCCGCCA'),('CAA-gain G1','CAACAACAGCCGCCA'),('CAA-gain G3–G8','CAA(CAG)k-1CAACAGCCGCCA'),('unclassified','other sequence')]
yy=Y1+H1-7.7
for i,(a,b) in enumerate(cfg):
    y_=yy-i*2.15; ax.add_patch(Rectangle((x+1.3,y_-0.55),0.9,1.1,fc=ORG,ec='none',zorder=4))
    T(x+3.0,y_,a,7.0,INK,'bold'); T(x+bw-1.0,y_,b,6.2,MUT,ha='right',family='DejaVu Sans Mono')
T(x+1.2,Y1+1.4,'Reads with an unresolved junction are not used',6.6,MUT,style='italic')
# arrows row 1
for i in range(3): arrow(X[i]+bw+0.25,Y1+H1/2,X[i+1]-0.25,Y1+H1/2)
# connector row1 -> row2
top2=Y2+H2; cy=top2+2.7
cxr=X[3]+bw/2
ax.plot([cxr,cxr],[Y1-0.2,cy],lw=1.6,color=MUT,zorder=5,solid_capstyle='butt')
ax.plot([cxr,X[0]+bw/2],[cy,cy],lw=1.6,color=MUT,zorder=5,solid_capstyle='butt')
arrow(X[0]+bw/2,cy,X[0]+bw/2,top2+0.2,MUT,1.6)
T(W/2,cy+1.3,'reads with a resolved structure: uninterrupted CAG count, configuration and CCG count',6.9,MUT,'normal','center',style='italic')
# ------------------------------------------------ row 2
# 5 Allele discovery
x=X[0]; box(x,Y2,bw,H2,'5','Allele discovery')
yt=Y2+H2
T(x+1.2,yt-6.0,'Modes of the CAG distribution, searched',7.4); T(x+1.2,yt-7.5,'separately within each configuration',7.4)
T(x+1.2,yt-9.4,'Score: reads within ±1 CAG (±2 above 35)',7.0,MUT); T(x+1.2,yt-10.9,'Called if ≥5% of resolved reads and ≥10 reads',7.0,MUT); T(x+1.2,yt-12.4,'At most two alleles',7.0,MUT)
hb=Y2+2.4
for i,v in enumerate([.12,.3,.6,1.0,.55,.25,.1]): ax.add_patch(Rectangle((x+1.8+i*1.05,hb),0.85,v*4.6,fc=BLUE,ec='none',zorder=4))
for i,v in enumerate([.08,.18,.35,.6,.9,1.0,.7,.3,.1]): ax.add_patch(Rectangle((x+10.6+i*1.05,hb),0.85,v*4.6,fc=BLUE,ec='none',zorder=4,alpha=0.55))
ax.plot([x+1.4,x+bw-1.0],[hb,hb],lw=0.8,color=MUT,zorder=5)
# 6 Stutter-aware sizing
x=X[1]; box(x,Y2,bw,H2,'6','Stutter-aware sizing (≥36 CAG)')
T(x+1.2,yt-6.0,'PCR stutter makes n−1 almost as abundant',7.4); T(x+1.2,yt-7.5,'as n: the allele is the longest count with',7.4); T(x+1.2,yt-9.0,'≥80% of the reads of the most frequent count',7.4)
bx0=x+3.0; base=Y2+3.3; sc=5.6
for i,(l,v,c) in enumerate([('n−2',.66,GRY),('n−1',1.0,GRY),('n',.96,BLUE),('n+1',.34,GRY)]):
    ax.add_patch(Rectangle((bx0+i*4.0,base),3.1,v*sc,fc=c,ec='none',zorder=4)); T(bx0+i*4.0+1.55,base-0.95,l,6.8,INK,w='bold' if l=='n' else 'normal',ha='center')
ax.plot([bx0-0.7,bx0+16.0],[base+0.8*sc]*2,ls=(0,(4,2)),lw=1.0,color=RED,zorder=6); T(bx0+16.3,base+0.8*sc,'80%',6.8,RED,'bold')
T(bx0+4.0+1.55,base+sc+0.8,'most frequent',6.4,MUT,ha='center'); T(bx0+8.0+1.55,base+0.96*sc+0.8,'called',6.8,BLUE,'bold','center')
# 7 Phasing
x=X[2]; box(x,Y2,bw,H2,'7','Phasing and fragment equivalent')
T(x+1.2,yt-6.0,'CCG length and configuration of each allele',7.4); T(x+1.2,yt-7.5,'= majority of its reads; chimeric reads ignored',7.4)
rows=[('Allele 1','17 CAG','canonical','CCG 7'),('Allele 2','39 CAG','canonical','CCG 7')]
for i,r in enumerate(rows):
    yy=yt-10.9-i*2.9
    ax.add_patch(FancyBboxPatch((x+1.3,yy-1.15),bw-2.6,2.3,boxstyle='round,pad=0,rounding_size=0.3',fc=PALE,ec=LINE,lw=0.7,zorder=3))
    T(x+1.9,yy,r[0],6.8,INK,'bold'); 
    for cc,(lab,colr,dx) in enumerate([(r[1],BLUE,6.0),(r[2],ORG,11.4),(r[3],GRN,17.0)]): T(x+dx+0.0,yy,lab,6.6,colr,'bold')
T(x+1.2,Y2+4.9,'Fragment-equivalent CAG = codons from the',7.0,MUT); T(x+1.2,Y2+3.4,'5′ anchor to the CCG-CCA motif − 4',7.0,MUT)
# 8 Output and flags
x=X[3]; box(x,Y2,bw,H2,'8','Output and flags')
T(x+1.2,yt-5.9,'summary.tsv · result.json · haplotype counts',7.2); T(x+1.2,yt-7.4,'HTML report with CAG, CAG×CCG and read plots',7.2)
hline(x+1.2,x+bw-1.2,yt-8.9)
fl=['Call differs from most frequent count (±1)','Alleles one repeat apart not separated','Marked allelic imbalance (ratio < 0.33)','Unclassified sequence: confirm independently']
for i,s in enumerate(fl):
    yy=yt-10.8-i*2.55; ax.add_patch(Polygon([[x+1.3,yy-0.55],[x+2.5,yy-0.55],[x+1.9,yy+0.65]],fc='#e0a21a',ec='none',zorder=4)); T(x+3.1,yy,s,6.8,INK)
# arrows row 2
for i in range(3): arrow(X[i]+bw+0.25,Y2+(H2)/2,X[i+1]-0.25,Y2+(H2)/2)
fig.savefig(OUT+'.pdf',metadata={'Creator':None,'Producer':None}); fig.savefig(OUT+'.png',dpi=170)
print(FONT)
