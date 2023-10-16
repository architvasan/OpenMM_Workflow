package require psfgen
package require autoionize
topology top_all36_prot.rtf
mol new 8gc_with_loops.pdb 
set all [atomselect top "protein"]
$all writepdb protein.pdb
pdbalias residue HIS HSE
pdbalias atom ILE CD1 CD
segment U {
    first ACE
    last CT2
    pdb protein.pdb
    }
coordpdb protein.pdb U
guesscoord
writepdb protein.gen.pdb
writepsf protein.gen.psf

package require solvate

solvate protein.gen.psf protein.gen.pdb -t 10 -o protein_wb

package require autoionize

autoionize -psf protein_wb.psf -pdb protein_wb.pdb -sc 0.15 -cation SOD -o apo
