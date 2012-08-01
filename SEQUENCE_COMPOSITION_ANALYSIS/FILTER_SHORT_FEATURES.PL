
# a following-ish command can be used to filter sequences shorter than 200bps.

awk '{ if( ($3-$2)>199 )print $0   }' wgEncodeUwDnaseH1hescHotspotsRep1.bed > wgEncodeUwDnaseH1hescHotspotsRep1_no_shorther_than_200bp.bed



