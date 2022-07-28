#include <stdio.h>
#include "hocdec.h"
#define IMPORT extern __declspec(dllimport)
IMPORT int nrnmpi_myid, nrn_nobanner_;

extern void _ca_reg();
extern void _cad_reg();
extern void _chanrhod_reg();
extern void _gabaa_reg();
extern void _kca_reg();
extern void _km_reg();
extern void _kv_reg();
extern void _na_reg();
extern void _na12_reg();
extern void _na16_reg();
extern void _ostim_reg();
extern void _xtra_reg();

void modl_reg(){
	//nrn_mswindll_stdio(stdin, stdout, stderr);
    if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
	fprintf(stderr, "Additional mechanisms from files\n");

fprintf(stderr," ca.mod");
fprintf(stderr," cad.mod");
fprintf(stderr," chanrhod.mod");
fprintf(stderr," gabaa.mod");
fprintf(stderr," kca.mod");
fprintf(stderr," km.mod");
fprintf(stderr," kv.mod");
fprintf(stderr," na.mod");
fprintf(stderr," na12.mod");
fprintf(stderr," na16.mod");
fprintf(stderr," ostim.mod");
fprintf(stderr," xtra.mod");
fprintf(stderr, "\n");
    }
_ca_reg();
_cad_reg();
_chanrhod_reg();
_gabaa_reg();
_kca_reg();
_km_reg();
_kv_reg();
_na_reg();
_na12_reg();
_na16_reg();
_ostim_reg();
_xtra_reg();
}
