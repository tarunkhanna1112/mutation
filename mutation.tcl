#  THIS CODE TRY TO REMOVE THE DEPENDENCE OF THE CODE ON THE INITIAL PDB FILE
proc builder { pdb pivot initial final k1 cresn flag inp} {

	set parm -1.0

	if  { $parm == "NA" } {
		set flag 1
	}

	package require math::linearalgebra

	# space variables

	set p(1) "   "
	set p(2) "  "
	set p(3) " "
	set p(4) ""
	set p(5) ""
	
	set p1(1) "    "
	set p1(2) "   "
	set p1(3) "  "
	set p1(4) " "
	set p1(5) ""

	set c(4) "    "
	set c(5) "   "
	set c(6) "  "
	set c(7) " "
	set c(8) ""

	set sat(1) "  "
	set sat(2) " "
	set sat(3) ""
	set sat(4) ""

	set ic(1) " "
	set ic(2) " "
	set ic(3) " "
	set ic(4) ""
	
	set satn(1) "   "
	set satn(2) "  "
	set satn(3) " "
	set satn(4) " "

	set f [open "$inp" "r"]
	set data [read $f]
	close $f

	set g [open "$pdb" "r"]
	set data1 [read $g]
	close $g

	set k 0

	while { [lindex $data $k] != $pivot } { 
		incr k
	}

	set xori1 [lindex $data [expr { $k + 4 }]]
	set yori1 [lindex $data [expr { $k + 5 }]]
	set zori1 [lindex $data [expr { $k + 6 }]]

	set k 0

	while { [lindex $data1 $k] != $pivot } {
		incr k
	}

	set cres [lindex $data1 [expr { $k - 1 }]]

	if { $initial == "NA" } {
		set initial $cres
	}

	set xori2 [lindex $data1 [expr { $k + 4 }]]
	set yori2 [lindex $data1 [expr { $k + 5 }]]
	set zori2 [lindex $data1 [expr { $k + 6 }]]

	set xshift [expr { $xori2 - $xori1 }]
	set yshift [expr { $yori2 - $yori1 }]
	set zshift [expr { $zori2 - $zori1 }]

	set h [open "int.pdb" "w"]

	if { $flag == 0 } {
		set k1 0
		while { [lindex $data $k1] != "C32" } {
			incr k1
		} 
			
		set xv [lindex $data [expr { $k1 + 4 }]]
		set xv [expr { $xv - $xori1 }]
		set yv [lindex $data [expr { $k1 + 5 }]]
		set yv [expr { $yv - $yori1}]
		set zv [lindex $data [expr { $k1 + 6 }]]
		set zv [expr { $zv - $zori1}]
	
		set v2 [list $xv $yv $zv]
		set nv2 [::math::linearalgebra::unitLengthVector $v2]

		set k1 0
		while { [lindex $data1 $k1] != "C32" } {
			incr k1
		} 

		set xv [lindex $data1 [expr { $k1 + 4 }]]
		set xv [expr { $xv - $xori2 }]
		set yv [lindex $data1 [expr { $k1 + 5 }]]
		set yv [expr { $yv - $yori2 }]
		set zv [lindex $data1 [expr { $k1 + 6 }]]
		set zv [expr { $zv - $zori2 }]
		
		set v1 [list $xv $yv $zv]
		set trans [::math::linearalgebra::sub $v2 $v1]
		set nv1 [::math::linearalgebra::unitLengthVector $v1]

		set angle [::math::linearalgebra::dotproduct $nv1 $nv2]
		set normal [::math::linearalgebra::crossproduct $nv2 $nv1]
		set normal [::math::linearalgebra::unitLengthVector $normal]

		set dum1 [open "dummy1" "w"]
		puts $dum1 "$normal"
		close $dum1

		set dum1 [open "dummy1" "r"]
		set du [read $dum1]
		close $dum1

		# ROTATION OF THE HEAD GROUP PART SPECIFIC TO THE FINAL LIPID ABOUT AN ARBITARY LINE PASSING THROUGH C31

		set a $xori1
		set b $yori1
		set cc $zori1
	
		set u [lindex $du 0]
		set v [lindex $du 1]
		set w [lindex $du 2]
		set u2 [expr { $u * $u }]
		set v2 [expr { $v * $v }]
		set w2 [expr { $w * $w }]
		set l [expr { $u2 + $v2 + $w2 }]

		set rl [expr { sqrt($l) }]
	
		#set theta [expr { ($angle * 3.14 * 0.0) / 180.0 }]
		
		set theta [expr { acos($angle) *  $parm }]
		#puts "[expr { 180 * $theta / 3.14 }]"

		set a11 [expr { ($u2 + (($v2 + $w2) * cos($theta))) / $l } ]
		set a12 [expr { ((($u*$v)*(1-cos($theta))) - ($w * $rl * sin($theta))) / $l }]
		set a13 [expr { ((($u*$w)*(1-cos($theta))) + ($v * $rl * sin($theta))) / $l }]
		set a21 [expr { ((($u*$v)*(1-cos($theta))) + ($w * $rl * sin($theta))) / $l }]
		set a22 [expr { ($v2 + (($u2 + $w2) * cos($theta))) / $l }]
		set a23 [expr { ((($v*$w)*(1-cos($theta))) - ($u * $rl * sin($theta))) / $l }]
		set a31 [expr { ((($u*$w)*(1-cos($theta))) - ($v * $rl * sin($theta))) / $l }]
		set a32 [expr { ((($v*$w)*(1-cos($theta))) + ($u * $rl * sin($theta))) / $l }]
		set a33 [expr { ($w2 + (($u2 + $v2) * cos($theta))) / $l }]
		set a41 [expr { (((($a * ($v2 + $w2)) - ($u * (($b*$v) + ($cc*$w)))) * (1 - cos($theta))) + (((($b*$w) - ($cc*$v))) * ($rl * sin($theta)))) / $l }]
		set a42 [expr { (((($b * ($u2 + $w2)) - ($v * (($a*$u) + ($cc*$w)))) * (1 - cos($theta))) + (((($cc*$u) - ($a*$w))) * ($rl * sin($theta)))) / $l }]
		set a43 [expr { (((($cc * ($u2 + $v2)) - ($w * (($a*$u) + ($b*$v)))) * (1 - cos($theta))) + (((($a*$v) - ($b*$u))) * ($rl * sin($theta)))) / $l }]

		set row1 [list $a11 $a12 $a13]
		set row2 [list $a21 $a22 $a23]
		set row3 [list $a31 $a32 $a33]

		set tr [list $row1 $row2 $row3]
		set tr1 [list $a41 $a42 $a43]
		#set tr [::math::linearalgebra::transpose $tr]

		# TRANSFORMATION 2 Z AXIS

		set tm2r1 [list [expr { cos($theta) }] [expr { -1 * sin($theta) }] 0.0] 
		set tm2r2 [list [expr { (sin($theta))}] [expr { cos($theta) }] 0.0]
		set tm2r3 [list 0.0 0.0 1.0]
		set tm2 [list $tm2r1 $tm2r2 $tm2r3]
		
				
	}

	set k 0
	while { $k < [llength $data1] } {
		if { [lindex $data1 $k] == "HETATM" } {
			set rn1 [lindex $data1 [expr { $k + 1 }]]
			set srn1 [string length $rn1]

			set at1 [lindex $data1 [expr { $k + 2 }]]
			set sat1 [string length $at1]

			set resn [lindex $data1 [expr { $k + 3 }]]
			set sresn [string length [lindex $data1 [expr { $k + 3 }]]]
			set sresn [expr { $sresn + 1 }]

			set cn [lindex $data1 [expr { $k + 4 }]]

			set an1 [lindex $data1 [expr { $k + 5 }]]
			set san1 [string length $an1]
			if { $san1 >= 4 } {
				set shift 1
				set an1 [string range [lindex $data1 [expr { $k + 4 }]] 1 4]
				set san1 4
				set cn a
			} else {
				set shift 0
			}

			if { $resn == "$cresn" && $rn1 >= $initial && $rn1 <= $final } {
				if { $flag == 0 } {
					set k1 0
				}
				while { [lindex $data $k1] != "$at1" } {
					puts "$at1 $k1"
					incr k1
				}

				set x1 [lindex $data [expr { $k1 + 4 - $shift}]]

				set y1 [lindex $data [expr { $k1 + 5 - $shift}]] 

				set z1 [lindex $data [expr { $k1 + 6 -$shift }]]
			
			} else {
				set x1 [lindex $data1 [expr { $k + 6 - $shift}]]
				set x1 [expr { $x1 - $xshift }]
				
				set y1 [lindex $data1 [expr { $k + 7 - $shift}]] 
				set y1 [expr { $y1 - $yshift }]

				set z1 [lindex $data1 [expr { $k + 8 -$shift }]]
				set z1 [expr { $z1 - $zshift }] 
			}
		
			if { $flag == 0 && $resn == $cresn && $rn1 < $initial} {

				set Rx [expr { $x1 - $xori1 }]
				set Ry [expr { $y1 - $yori1 }]
				set Rz [expr { $z1 - $zori1 }]
				set R [expr { ($Rx*$Rx) + ($Ry*$Ry) + ($Rz*$Rz) }]
				set R [expr { sqrt($R) }]
				set tori [list $xori1 $yori1 $zori1]
				set tcoord [list $x1 $y1 $z1]

				set tcoord [::math::linearalgebra::matmul $tr $tcoord]
				set tcoord [::math::linearalgebra::add $tcoord $tr1]
				set init $tcoord
	
				set dum [open "dummy" "w"]

				puts $dum "$tcoord"

				close $dum

				set dum [open "dummy" "r"]
				set dvar [read $dum]
				close $dum

				set x1 [format "%.3f" [lindex $dvar 0]]
				set sx1 [string length $x1]
	
				set y1 [format "%.3f" [lindex $dvar 1]]
				set sy1 [string length $y1]
				
				set z1 [format "%.3f" [lindex $dvar 2]]
				set sz1 [string length $z1]

				
			}							
			set x1 [format "%.3f" $x1]
			set sx1 [string length $x1]
			set y1 [format "%.3f" $y1]
			set sy1 [string length $y1]
			set z1 [format "%.3f" $z1]
			set sz1 [string length $z1]
			puts $h "HETATM$p1($srn1)$rn1 $ic($sat1)$at1$satn($sat1)$resn $p($sresn)a$p($san1)$an1    $c($sx1)$x1$c($sy1)$y1$c($sz1)$z1  1.00  0.00           [lindex $data1 [expr { $k + 11 - $shift }]]"
		}
		incr k
	}
	puts $h "END"
	close $h
	return $k1
}

set OL 51
set PA 46
set MY 40
set LA 34

if { [lindex $::argv 5] == "OL" } {
	set P $OL 
}
if { [lindex $::argv 5] == "PA" } {
	set P $PA 
}
if { [lindex $::argv 5] == "MY" } {
	set P $MY 
}
if { [lindex $::argv 5] == "LA" } {
	set P $LA 
}

set f [open "database" "r"]
set data [read $f]
close $f

set t1 [lindex $::argv 5]
set k 0
while { [lindex $data 1 $k] != $t1 } {
	incr k
}
set inp1 [lindex $data 0 $k]
set t2 [lindex $::argv 7]
set k 0
while { [lindex $data 1 $k] != $t2 } {
	incr k
}
set inp2 [lindex $data 0 $k]

if { $inp1 == $inp2 } {
	puts "				**** NEW PDB FORMED OUT OF $inp1 ****"
} else { 
	puts "				**** NEW PDB FORMED OUT OF $inp1 AND $inp2 ****"
}

# HEAD
set d0 [builder [lindex $::argv 1] C31 NA 140 0 [lindex $::argv 3] 0 $inp1]
# TAIL 1
set d [builder int.pdb C11 1 $P 0 $t1 1 $inp1]
# TAIL 2
set d2 [builder int.pdb C21 $P 140 $d $t2 1 $inp2]

#exec mkdir -p new_pdb

exec mv int.pdb new_[lindex $::argv 1]

#exec cp [lindex $::argv 1]_new.pdb ./new_pdb





