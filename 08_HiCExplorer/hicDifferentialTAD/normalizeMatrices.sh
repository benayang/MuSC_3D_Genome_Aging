# hicNormalize --matrices aged.merged_40000.cool young.merged_40000.cool \
# --normalize smallest \
# --outFileName aged.merged_40000_minNorm.cool young.merged_40000_minNorm.cool

hicConvertFormat -m aged.merged_40000_minNorm.cool \
--inputFormat cool --outputFormat cool \
-o aged.merged_40000_minNorm_KR.cool \
--correction_name KR

hicConvertFormat -m young.merged_40000_minNorm.cool \
--inputFormat cool --outputFormat cool \
-o young.merged_40000_minNorm_KR.cool \
--correction_name KR