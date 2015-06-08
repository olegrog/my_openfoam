#!/bin/bash

. $WM_PROJECT_DIR/bin/tools/RunFunctions

cases=$(ls | grep '^[0-9]' | sort -g)

norm() {
    echo "print $1/$2" | python
}

for U in $cases; do
    echo "U = $U"
    (
    printf '#%11s %12s %11s %12s %11s %11s %11s %11s %11s %11s %11s %12s %11s %11s\n' \
        Kn Pxy M Qx Qy Pxx Pyy Pzz tau P DU0 DT0 p0 T_B0
    for Kn in $(ls $U | sort -g); do
        last=$(grep ^Time $U/$Kn/log.couetteFlowFoam | tail -1 | awk '{ print $3 }')
        DU=$(grep -A4 top $U/$Kn/$last/wallDU0 | grep value | sed 's/.* //;s/;//')
        DT=$(grep -A4 top $U/$Kn/$last/wallDT0 | grep value | sed 's/.* //;s/;//')
        TB=$(grep -A4 top $U/$Kn/$last/T0 | grep value | sed 's/.* //;s/;//')
        p0=$(grep uniform $U/$Kn/$last/press0 | sed 's/.* //;s/;//')
        M=$(norm $(tail $U/$Kn/log.couetteFlowFoam | grep 'M =' | awk '{print $3}') $U)
        tau=$(norm $(tail $U/$Kn/log.couetteFlowFoam | grep 'tau =' | awk '{print $3}') $U)
        P=$(norm $(tail $U/$Kn/log.couetteFlowFoam | grep 'P =' | awk '{print $3}') $U)
        Pxy=$(norm $(tail $U/$Kn/log.couetteFlowFoam | grep 'P_ij =' | awk '{print $4}') $U)
        Pxx=$(norm $(tail $U/$Kn/log.couetteFlowFoam | grep 'P_ij =' | awk '{print $3}' | sed 's/(//') $U)
        Pyy=$(norm $(tail $U/$Kn/log.couetteFlowFoam | grep 'P_ij =' | awk '{print $6}') $U)
        Pzz=$(norm $(tail $U/$Kn/log.couetteFlowFoam | grep 'P_ij =' | awk '{print $8}' | sed 's/)//') $U)
        Qx=$(norm $(tail $U/$Kn/log.couetteFlowFoam | grep 'q_i =' | awk '{print $3}' | sed 's/(//') $U)
        Qy=$(norm $(tail $U/$Kn/log.couetteFlowFoam | grep 'q_i =' | awk '{print $4}') $U)
        printf "%.6e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e\n" $Kn $Pxy $M $Qx $Qy $Pxx $Pyy $Pzz $tau $P $DU $DT $p0 $TB
    done
    ) > $1-$(printf "%.1f" $U).txt
done
