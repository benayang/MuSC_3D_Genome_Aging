prefix=/nas/homes/benyang/HiC

export _JAVA_OPTIONS="-Xms256M -Xmx50g"

juicerJar=$prefix/juicer_tools_1.22.01.jar

# https://groups.google.com/g/3d-genomics/c/LSY9NSfsyX8

$juicerJar apa -u $prefix/02_HIC/aged.merged/aged.merged.hic \
$prefix/05_Arrowhead_HICCUPS/aged_merged_hiccups/merged_loops.bedpe \
$prefix/05_Arrowhead_HICCUPS/aged_merged_APA \
--threads 1

$juicerJar apa -u $prefix/02_HIC/young.merged/young.merged.hic \
$prefix/05_Arrowhead_HICCUPS/young_merged_hiccups/merged_loops.bedpe \
$prefix/05_Arrowhead_HICCUPS/young_merged_APA \
--threads 1

## Run APA on differential loops

$juicerJar apa -u $prefix/02_HIC/young.merged/young.merged.hic \
$prefix/05_Arrowhead_HICCUPS/hiccups_diff/young_differential_loops1.bedpe \
$prefix/05_Arrowhead_HICCUPS/hiccups_diff_APA/young_diff_young_HiC \
--threads 1

$juicerJar apa -u $prefix/02_HIC/young.merged/young.merged.hic \
$prefix/05_Arrowhead_HICCUPS/hiccups_diff/aged_differential_loops2.bedpe \
$prefix/05_Arrowhead_HICCUPS/hiccups_diff_APA/aged_diff_young_HiC \
--threads 1

$juicerJar apa -u $prefix/02_HIC/aged.merged/aged.merged.hic \
$prefix/05_Arrowhead_HICCUPS/hiccups_diff/young_differential_loops1.bedpe \
$prefix/05_Arrowhead_HICCUPS/hiccups_diff_APA/young_diff_aged_HiC \
--threads 1

$juicerJar apa -u $prefix/02_HIC/aged.merged/aged.merged.hic \
$prefix/05_Arrowhead_HICCUPS/hiccups_diff/aged_differential_loops2.bedpe \
$prefix/05_Arrowhead_HICCUPS/hiccups_diff_APA/aged_diff_aged_HiC \
--threads 1