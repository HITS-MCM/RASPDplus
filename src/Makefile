

BIN=pdb2chaininfo pdb2noh pdb2pdb protein-res-cm-cal CM-PL print_volm_1M print_volm_customized protein_param_calc paste final_select_all final_select maxD_calculator pdbtopdb windex_single lp1 lp2 lp3 lp5 lp6 lp7 lp8 hbd_hba lipi paste_ml_format paste_ml_format_4smiles print_volm_multiple pdbarrange mass param_chek_nonzero

all: $(BIN) connectmol

connectmol: 
	cd CONNECTMOLv2_t2 && $(MAKE) 
	cp CONNECTMOLv2_t2/Connect2.0.exe ../bin/	

pdb2chaininfo: pdb2chaininfo.c
	g++ pdb2chaininfo.c -o ../bin/pdb2chaininfo.exe

pdb2noh: pdb2noh.c
	g++ pdb2noh.c -o ../bin/pdb2noh.exe

pdb2pdb: pdb2pdb.c
	g++ pdb2pdb.c -o ../bin/pdb2pdb.exe

protein-res-cm-cal: protein-res-cm-cal.c
	g++ protein-res-cm-cal.c -o ../bin/protein-res-cm-cal.exe

CM-PL: CM-PL.c
	g++ CM-PL.c -o ../bin/CM-PL.exe

print_volm_1M: print_volm_1M.c
	g++ print_volm_1M.c -o ../bin/print_volm_1M.exe

print_volm_customized: print_volm_customized.c
	g++ print_volm_customized.c -o ../bin/print_volm_customized.exe

protein_param_calc: protein_param_calc.cpp
	g++ protein_param_calc.cpp -o ../bin/protein_param_calc.exe

paste: paste.c
	g++ paste.c -o ../bin/paste.exe

final_select_all: final_select_all.c
	g++ final_select_all.c -o ../bin/final_select_all.exe

final_select: final_select.c
	g++ final_select.c -o ../bin/final_select.exe

maxD_calculator: maxD_calculator.cpp
	g++ maxD_calculator.cpp -o ../bin/maxD_calculator.exe

pdbtopdb: pdbtopdb.c
	g++ pdbtopdb.c -o ../bin/pdbtopdb.exe

windex_single: windex_single.cpp
	g++ windex_single.cpp -o ../bin/windex_single.exe

lp1: lp1.c
	g++ lp1.c -o ../bin/lp1.exe

lp2: lp2.c
	g++ lp2.c -o ../bin/lp2.exe

lp3: lp3.c
	g++ lp3.c -o ../bin/lp3.exe

lp5: lp5.c
	g++ lp5.c -o ../bin/lp5.exe

lp6: lp6.c
	g++ lp6.c -o ../bin/lp6.exe

lp7: lp7.c
	g++ lp7.c -o ../bin/lp7.exe

lp8: lp8.c
	g++ lp8.c -o ../bin/lp8.exe

hbd_hba: hbd_hba.c
	g++ hbd_hba.c -o ../bin/hbd_hba.exe

lipi: lipi.c
	g++ lipi.c -o ../bin/lipi.exe

paste_ml_format: paste_ml_format.c
	g++ paste_ml_format.c -o ../bin/paste_ml_format.exe

paste_ml_format_4smiles: paste_ml_format_4smiles.c
	g++ paste_ml_format_4smiles.c -o ../bin/paste_ml_format_4smiles.exe

print_volm_multiple: print_volm_multiple.c
	g++ print_volm_multiple.c -o ../bin/print_volm_multiple.exe

pdbarrange: pdbarrange.c
	g++ pdbarrange.c -o ../bin/pdbarrange.exe

mass: mass.c
	g++ mass.c -o ../bin/mass.exe

param_chek_nonzero: param_chek_nonzero.c
	g++ param_chek_nonzero.c -o ../bin/param_chek_nonzero.exe
