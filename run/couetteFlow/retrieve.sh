#!/bin/bash

. $WM_PROJECT_DIR/bin/tools/RunFunctions

cases=$(ls | grep '^[0-9]' | sort -g)
if test $# -gt 0; then
    cases=$@
fi

norm() {
    echo "print $1/$2" | python
}

for U in $cases; do
    echo "U = $U"
    (
    printf '#%11s %13s %12s %13s %12s %12s %12s %12s %12s %12s %13s %12s %12s\n' \
        Kn Pxy M Qx Qy Pxx Pyy Pzz T DU0 DT0 p0 T_B0
    for Kn in $(ls $U | sort -g); do
        last=$(grep ^Time $U/$Kn/log.couetteFlowFoam | tail -1 | awk '{ print $3 }')
        DU=$(grep -A4 top $U/$Kn/$last/wallDU0 | grep value | sed 's/.* //;s/;//')
        DT=$(grep -A4 top $U/$Kn/$last/wallDT0 | grep value | sed 's/.* //;s/;//')
        TB=$(grep -A4 top $U/$Kn/$last/T0 | grep value | sed 's/.* //;s/;//')
        p0=$(grep uniform $U/$Kn/$last/p0 | sed 's/.* //;s/;//')
        M=$(norm $(tail $U/$Kn/log.couetteFlowFoam | grep 'M =' | awk '{print $3}') $U)
        T=$(norm $(tail $U/$Kn/log.couetteFlowFoam | grep 'T =' | awk '{print $3}') $U)
        Pxx=$(norm $(tail $U/$Kn/log.couetteFlowFoam | grep 'p_ij =' | awk '{print $3}' | sed 's/(//') $U)
        Pxy=$(norm $(tail $U/$Kn/log.couetteFlowFoam | grep 'p_ij =' | awk '{print $4}') "($U*0.5)")
        Pyy=$(norm $(tail $U/$Kn/log.couetteFlowFoam | grep 'p_ij =' | awk '{print $6}') $U)
        Pzz=$(norm $(tail $U/$Kn/log.couetteFlowFoam | grep 'p_ij =' | awk '{print $8}' | sed 's/)//') $U)
        Qx=$(norm $(tail $U/$Kn/log.couetteFlowFoam | grep 'q_i =' | awk '{print $3}' | sed 's/(//') $U)
        Qy=$(norm $(tail $U/$Kn/log.couetteFlowFoam | grep 'q_i =' | awk '{print $4}') $U)
        printf "%.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e\n" $Kn $Pxy $M $Qx $Qy $Pxx $Pyy $Pzz $T $DU $DT $p0 $TB
    done
    ) > ns-$(printf "%.1f" $U).txt
done
