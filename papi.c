#include "XSbench_header.h"

void counter_init( int * eventset )
{
	int names[NUM_PAPI_EVENTS] = 
	{
		PAPI_TOT_INS,
		PAPI_INT_INS,
		PAPI_FP_INS,
		PAPI_LD_INS,
		PAPI_SR_INS,
		PAPI_BR_INS,
		PAPI_VEC_INS,
		PAPI_RES_STL,
		PAPI_FP_STAL,
		PAPI_TOT_CYC,
		PAPI_LST_INS,
		PAPI_SYC_INS,
		PAPI_L1_DCH,
		PAPI_L2_DCH,
		PAPI_L1_DCA,
		PAPI_L2_DCA,
		PAPI_L3_DCA,
		PAPI_L1_DCR,
		PAPI_L2_DCR,
		PAPI_L3_DCR,
		PAPI_L1_DCW,
		PAPI_L2_DCW,
		PAPI_L3_DCW,
		PAPI_L1_ICH,
		PAPI_L2_ICH,
		PAPI_L3_ICH,
		PAPI_L1_ICA,
		PAPI_L2_ICA,
		PAPI_L3_ICA,
		PAPI_L1_ICR,
		PAPI_L2_ICR,
		PAPI_L3_ICR,
		PAPI_L1_ICW,
		PAPI_L2_ICW,
		PAPI_L3_ICW,
		PAPI_L1_TCH,
		PAPI_L2_TCH,
		PAPI_L3_TCH,
		PAPI_L1_TCA,
		PAPI_L2_TCA,
		PAPI_L3_TCA,
		PAPI_L1_TCR,
		PAPI_L2_TCR,
		PAPI_L3_TCR,
		PAPI_L1_TCW,
		PAPI_L2_TCW,
		PAPI_L3_TCW,
		PAPI_FML_INS,
		PAPI_FAD_INS,
		PAPI_FDV_INS,
		PAPI_FSQ_INS,
		PAPI_FNV_INS,
		PAPI_FP_OPS,
		PAPI_SP_OPS,
		PAPI_DP_OPS,
		PAPI_VEC_SP,
		PAPI_VEC_DP,
		PAPI_REF_CYC,
	};
	
	PAPI_library_init(83886080);

	PAPI_create_eventset(eventset);

	PAPI_add_events( *eventset, names, NUM_PAPI_EVENTS );

	PAPI_start(*eventset);
}

void counter_stop( int * eventset )
{
	long_long values[NUM_PAPI_EVENTS] = {0};

	PAPI_stop(*eventset, values);

	FILE * out = fopen("counters.txt", "w");

	for( int i = 0; i < NUM_PAPI_EVENTS; i++ )
	{
		fprintf(out,"%lld\n", values[i]);
	}
	
	fclose(out);
}
