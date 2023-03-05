import fanc
import fanc.plotting as fancplot
import numpy as np
import matplotlib.pyplot as plt

aged=fanc.load('/nas/homes/benyang/HiC/02_HIC/aged.merged/aged.merged.hic@5kb@KR')
young=fanc.load('/nas/homes/benyang/HiC/02_HIC/young.merged/young.merged.hic@5kb@KR')

aged_TAD=fancplot.GenePlot("/nas/homes/benyang/HiC/08_HiCExplorer/aged.merged/40kb/aged.merged_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_domains.bed",aspect=0.25)
young_TAD=fancplot.GenePlot("/nas/homes/benyang/HiC/08_HiCExplorer/young.merged/40kb/young.merged_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_domains.bed",aspect=0.25)

# young_AB=fancplot.GenePlot("/nas/homes/benyang/HiC/04_FANC/without_KR_normalization/young.ab_100kb.bed",aspect=0.25)
# aged_AB=fancplot.GenePlot("/nas/homes/benyang/HiC/04_FANC/without_KR_normalization/aged.ab_100kb.bed",aspect=0.25)

genes=fancplot.GenePlot("/nas/homes/benyang/HiC/get_tss/genebodies.gencode.vM25.basic.annotation.filtered.uniq.knownGenes.bed",aspect=0.4)

def saveData(gene):
    young_line_xd = plt.gca().get_lines()[0].get_xdata()
    young_line_yd = plt.gca().get_lines()[0].get_ydata()
    aged_line_xd = plt.gca().get_lines()[1].get_xdata()
    aged_line_yd = plt.gca().get_lines()[1].get_ydata()
    young_data = np.array([young_line_xd, young_line_yd])
    aged_data = np.array([aged_line_xd, aged_line_yd])
    np.savetxt(gene+"_young_data.txt", young_data, delimiter="\t")
    np.savetxt(gene+"_aged_data.txt", aged_data, delimiter="\t")

# mef2a_plot4C=fancplot.HicSlicePlot([young,aged],"chr7:67231163-67372858",colors=['#0072B5FF','#BC3C29FF'])
# # mef2a_plot4C=fancplot.HicSlicePlot([young,aged],"chr7:67226163-67236163",colors=['#0072B5FF','#BC3C29FF'])
# gfig=fancplot.GenomicFigure([mef2a_plot4C])
# fig,ax=gfig.plot("chr7:66653746-68376854")
# fig.savefig("virtual_4c_plot_Mef2a.png")
# saveData("Mef2a")

# mtor_plot4C=fancplot.HicSlicePlot([young,aged],"chr4:148448582-148557685",colors=['#0072B5FF','#BC3C29FF'])
# # mtor_plot4C=fancplot.HicSlicePlot([young,aged],"chr4:148443582-148453582",colors=['#0072B5FF','#BC3C29FF'])
# gfig=fancplot.GenomicFigure([mtor_plot4C])
# fig,ax=gfig.plot("chr4:147824771-149729874")
# fig.savefig("virtual_4c_plot_Mtor.png")
# saveData("Mtor")

# foxo3_plot4C=fancplot.HicSlicePlot([young,aged],"chr10:42185786-42276742",colors=['blue','red'])
# # foxo3_plot4C=fancplot.HicSlicePlot([young,aged],"chr10:42180786-42190786",colors=['blue','red'])
# gfig=fancplot.GenomicFigure([foxo3_plot4C])
# fig,ax=gfig.plot("chr10:40991533-42769626")
# fig.savefig("virtual_4c_plot_Foxo3.png")
# saveData("Foxo3")

# mxi1_plot4C=fancplot.HicSlicePlot([young,aged],"chr19:53310506-53375810",colors=['blue','red'])
# # mxi1_plot4C=fancplot.HicSlicePlot([young,aged],"chr19:53305506-53315506",colors=['blue','red'])
# gfig=fancplot.GenomicFigure([mxi1_plot4C])
# fig,ax=gfig.plot("chr19:50127905-54067414")
# fig.savefig("virtual_4c_plot_Mxi1.png")
# saveData("Mxi1")

# sesn2_plot4C=fancplot.HicSlicePlot([young,aged],"chr4:132492804-132510456",colors=['blue','red'])
# # sesn2_plot4C=fancplot.HicSlicePlot([young,aged],"chr4:132487804-132497804",colors=['blue','red'])
# gfig=fancplot.GenomicFigure([sesn2_plot4C])
# fig,ax=gfig.plot("chr4:131831552-133564320")
# fig.savefig("virtual_4c_plot_Sesn2.png")
# saveData("Sesn2")

# pax7_plot4C=fancplot.HicSlicePlot([young,aged],"chr4:139737059-139833026",colors=['blue','red'])
# gfig=fancplot.GenomicFigure([pax7_plot4C])
# fig,ax=gfig.plot("chr4:138868292-140531190")
# fig.savefig("virtual_4c_plot_Pax7.png")
# saveData("Pax7")

# prdm16_plot4C=fancplot.HicSlicePlot([young,aged],"chr4:154314122-154638873",colors=['blue','red'])
# gfig=fancplot.GenomicFigure([prdm16_plot4C])
# fig,ax=gfig.plot("chr4:153865383-155344873")
# fig.savefig("virtual_4c_plot_Prdm16.png")
# saveData("Prdm16")

# sesn1_plot4C=fancplot.HicSlicePlot([young,aged],"chr10:41808574-41910436",colors=['blue','red'])
# gfig=fancplot.GenomicFigure([sesn1_plot4C,young_TAD,aged_TAD,genes])
# fig,ax=gfig.plot("chr10:41000000-42500000")
# fig.savefig("virtual_4c_plot_Sesn1.png")

# ddit3_plot4C=fancplot.HicSlicePlot([young,aged],"chr10:127288793-127298288",colors=['#0072B5FF','#BC3C29FF'])
# gfig=fancplot.GenomicFigure([ddit3_plot4C])
# fig,ax=gfig.plot("chr10:126686384-128047303")
# fig.savefig("virtual_4c_plot_Ddit3.png")
# saveData("Ddit3")

# kdm5b_plot4C=fancplot.HicSlicePlot([young,aged],"chr1:134560178-134632878",colors=['#0072B5FF','#BC3C29FF'])
# gfig=fancplot.GenomicFigure([kdm5b_plot4C])
# fig,ax=gfig.plot("chr1:134113889-135443577")
# fig.savefig("virtual_4c_plot_Kdm5b.png")
# saveData("Kdm5b")

# cdkn1a_plot4C=fancplot.HicSlicePlot([young,aged],"chr17:29090979-29100722",colors=['#0072B5FF','#BC3C29FF'])
# gfig=fancplot.GenomicFigure([cdkn1a_plot4C,young_TAD,aged_TAD,genes])
# fig,ax=gfig.plot("chr17:28903012-29568771")
# fig.savefig("virtual_4c_plot_Cdkn1a.png")
# saveData("Cdkn1a")

# spry2_plot4C=fancplot.HicSlicePlot([young,aged],"chr14:105891947-105896819",colors=['#0072B5FF','#BC3C29FF'])
# gfig=fancplot.GenomicFigure([spry2_plot4C])
# fig,ax=gfig.plot("chr14:105201592-106403934")
# fig.savefig("virtual_4c_plot_Spry2.png")
# saveData("Spry2")

# stat1_plot4C=fancplot.HicSlicePlot([young,aged],"chr1:52119438-52161865",colors=['#0072B5FF','#BC3C29FF'])
# gfig=fancplot.GenomicFigure([stat1_plot4C])
# fig,ax=gfig.plot("chr1:51699207-52716926")
# fig.savefig("virtual_4c_plot_Stat1.png")
# saveData("Stat1")

# kcnn2_plot4C=fancplot.HicSlicePlot([young,aged],"chr18:45268860-45685883",colors=['#0072B5FF','#BC3C29FF'])
# gfig=fancplot.GenomicFigure([kcnn2_plot4C])
# fig,ax=gfig.plot("chr18:44560504-46449466")
# fig.savefig("virtual_4c_plot_Kcnn2.png")
# saveData("Kcnn2")

# pou3f1_plot4C=fancplot.HicSlicePlot([young,aged],"chr4:124657646-124660655",colors=['#0072B5FF','#BC3C29FF'])
# gfig=fancplot.GenomicFigure([pou3f1_plot4C])
# fig,ax=gfig.plot("chr4:124611663-124713519")
# fig.savefig("virtual_4c_plot_Pou3f1.png")
# saveData("Pou3f1")

mpc2_plot4C=fancplot.HicSlicePlot([young,aged],"chr1:165461208-165481214",colors=['#0072B5FF','#BC3C29FF'])
gfig=fancplot.GenomicFigure([mpc2_plot4C])
fig,ax=gfig.plot("chr1:165216530-165725389")
fig.savefig("virtual_4c_plot_Mpc2.png")
saveData("Mpc2")

mpc1_plot4C=fancplot.HicSlicePlot([young,aged],"chr17:8283813-8297661",colors=['#0072B5FF','#BC3C29FF'])
gfig=fancplot.GenomicFigure([mpc1_plot4C])
fig,ax=gfig.plot("chr17:8166577-8413850")
fig.savefig("virtual_4c_plot_mpc1.png")
saveData("Mpc1")
