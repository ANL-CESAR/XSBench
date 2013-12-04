#include "XSbench_header.h"

void counter_init( int *eventset, int *num_papi_events )
{
	char error_str[PAPI_MAX_STR_LEN];
	int stat;
	
	/////////////////////////////////////////////////////////////////////////
	//                        PAPI EVENT SELECTION
	/////////////////////////////////////////////////////////////////////////
	// User can comment/uncomment blocks as they see fit within this seciton

	// Some Standard Events
	//int events[] = {PAPI_TOT_INS,PAPI_LD_INS,PAPI_FP_INS};
	
	// Bandwidth Used
	// ((PAPI_Lx_TCM * Lx_linesize) / PAPI_TOT_CYC) * Clock(MHz)
	//int events[] = {PAPI_L3_TCM, PAPI_TOT_CYC};

	// L3 Total Cache Miss Ratio
	// PAPI_L3_TCM / PAPI_L3_TCA
	// (On Xeon dual octo -  65%, not dependent on # of threads)
	//int events[] = {PAPI_L3_TCM, PAPI_L3_TCA};
	
	// % Cycles with no instruction use
	// PAPI_STL_ICY / PAPI_TOT_CYC
	//int events[] = { PAPI_STL_ICY, PAPI_TOT_CYC };
	
	// % Branch instructions Mispredicted
	// PAPI_BR_MSP / PAPI_BR_CN
	//int events[] = { PAPI_BR_MSP, PAPI_BR_CN, PAPI_BR_PRC };
	 
	// TLB Misses
	//int events[] = { PAPI_TLB_DM };

	// MFlops
	// (PAPI_FP_INS/PAPI_TOT_CYC) * Clock(MHz)
	//int events[] = { PAPI_FP_INS, PAPI_TOT_CYC };

	
	// TLB misses (Using native counters)
	int events[2];
	int EventCode;
	char * event1 = "perf::DTLB-LOADS";
	char * event2 = "perf::DTLB-LOAD-MISSES";
	PAPI_event_name_to_code( event1, &EventCode );
	events[0] = EventCode;	
	PAPI_event_name_to_code( event2, &EventCode );
	events[1] = EventCode;	

	/*	
	// Stalled Cycles, front v back (Using native counters)
	int events[3];
	int EventCode;
	char * event1 = "perf::STALLED-CYCLES-FRONTEND";
	char * event2 = "perf::STALLED-CYCLES-BACKEND";
	char * event3 = "perf::PERF_COUNT_HW_CPU_CYCLES";
	PAPI_event_name_to_code( event1, &EventCode );
	events[0] = EventCode;	
	PAPI_event_name_to_code( event2, &EventCode );
	events[1] = EventCode;	
	PAPI_event_name_to_code( event3, &EventCode );
	events[2] = EventCode;	
	*/	
	/*
	// LLC Cache Misses (Using native counters)
	int events[2];
	int EventCode;
	char * event1 = "ix86arch::LLC_REFERENCES";
	char * event2 = "ix86arch::LLC_MISSES";
	PAPI_event_name_to_code( event1, &EventCode );
	events[0] = EventCode;	
	PAPI_event_name_to_code( event2, &EventCode );
	events[1] = EventCode;	
	*/

	/*
	// Node Prefetch Misses (Using native counters)
	int events[1];
	int EventCode;
	//char * event1 = "perf::NODE-PREFETCHES";
	//char * event2 = "perf::NODE-PREFETCH-MISSES";
	char * event1 = "perf::NODE-PREFETCHES";
	char * event2 = "perf::NODE-LOAD-MISSES:COUNT";
	//PAPI_event_name_to_code( event1, &EventCode );
	//events[0] = EventCode;	
	PAPI_event_name_to_code( event2, &EventCode );
	events[0] = EventCode;	
	*/

	/*
	// CPU Stalls Due to lack of Load Buffers (Using native counters)
	int events[2];
	int EventCode;
	char * event1 = "RESOURCE_STALLS:LB";
	char * event2 = "perf::PERF_COUNT_HW_CPU_CYCLES";
	PAPI_event_name_to_code( event1, &EventCode );
	events[0] = EventCode;	
	PAPI_event_name_to_code( event2, &EventCode );
	events[1] = EventCode;	
	*/	
	/*
	// CPU Stalls Due to ANY Resource (Using native counters)
	int events[2];
	int EventCode;
	char * event1 = "RESOURCE_STALLS:ANY";
	char * event2 = "perf::PERF_COUNT_HW_CPU_CYCLES";
	PAPI_event_name_to_code( event1, &EventCode );
	events[0] = EventCode;	
	PAPI_event_name_to_code( event2, &EventCode );
	events[1] = EventCode;	
	*/	

	/*	
	// CPU Stalls at Reservation Station (Using native counters)
	int events[2];
	int EventCode;
	char * event1 = "RESOURCE_STALLS:RS";
	char * event2 = "perf::PERF_COUNT_HW_CPU_CYCLES";
	PAPI_event_name_to_code( event1, &EventCode );
	events[0] = EventCode;	
	PAPI_event_name_to_code( event2, &EventCode );
	events[1] = EventCode;	
	*/

	// CPU Stall Reason Breakdown (Using native counters)
	/*
	int events[4];
	int EventCode;
	// Set 1
	char * event1 = "RESOURCE_STALLS:ANY";
	char * event2 = "RESOURCE_STALLS:LB";
	char * event3 = "RESOURCE_STALLS:RS";
	char * event4 = "RESOURCE_STALLS:SB";
	// Set 1
	// Set 2
	char * event1 = "RESOURCE_STALLS:ANY";
	char * event2 = "RESOURCE_STALLS:ROB";
	char * event3 = "RESOURCE_STALLS:MEM_RS";
	char * event4 = "RESOURCE_STALLS2:ALL_FL_EMPTY";
	// Set 2
	// Set 3
	char * event1 = "RESOURCE_STALLS:ANY";
	char * event2 = "RESOURCE_STALLS2:ALL_PRF_CONTROL";
	char * event3 = "RESOURCE_STALLS2:ANY_PRF_CONTROL"; // duplicate
	char * event4 = "RESOURCE_STALLS2:OOO_RSRC";
	// Set 3

	// Events that don't need to be counted
	// Don't bother measuring these
	//char * event1 = "RESOURCE_STALLS:FCSW"; // Always 0, don't measure
	//char * event1 = "RESOURCE_STALLS:MXCSR"; // Always 0, don't measure
	//char * event3 = "RESOURCE_STALLS2:BOB_FULL"; // Always trivial
	//char * event3 = "RESOURCE_STALLS2:ANY_PRF_CONTROL"; // duplicate
	
	PAPI_event_name_to_code( event1, &EventCode );
	events[0] = EventCode;	
	PAPI_event_name_to_code( event2, &EventCode );
	events[1] = EventCode;	
	PAPI_event_name_to_code( event3, &EventCode );
	events[2] = EventCode;	
	PAPI_event_name_to_code( event4, &EventCode );
	events[3] = EventCode;	
	*/

	/////////////////////////////////////////////////////////////////////////
	//                        PAPI EVENT LOADING
	/////////////////////////////////////////////////////////////////////////
	// Users should not need to alter anything within this section

	int thread = omp_get_thread_num();
	if( thread == 0 )
		printf("Initializing PAPI counters...\n");

	*num_papi_events = sizeof(events) / sizeof(int);

	if ((stat = PAPI_thread_init((long unsigned int (*)(void)) omp_get_thread_num)) != PAPI_OK){
		PAPI_perror("PAPI_thread_init");
		exit(1);
	}

	if ( (stat= PAPI_create_eventset(eventset)) != PAPI_OK){
		PAPI_perror("PAPI_create_eventset");
		exit(1);
	}

	for( int i = 0; i < *num_papi_events; i++ ){
		if ((stat=PAPI_add_event(*eventset,events[i])) != PAPI_OK){
			PAPI_perror("PAPI_add_event");
			exit(1);
		}
	}

	if ((stat=PAPI_start(*eventset)) != PAPI_OK){
		PAPI_perror("PAPI_start");
		exit(1);
	}
}

// Stops the papi counters and prints results
void counter_stop( int * eventset, int num_papi_events )
{
	int * events = malloc(num_papi_events * sizeof(int));
	int n = num_papi_events;
	PAPI_list_events( *eventset, events, &n );
	PAPI_event_info_t info;

	long_long * values = malloc( num_papi_events * sizeof(long_long));
	PAPI_stop(*eventset, values);
	int thread = omp_get_thread_num();

#pragma omp critical (papi)
	{
		printf("Thread %d\n", thread);
		for( int i = 0; i < num_papi_events; i++ )
		{
			PAPI_get_event_info(events[i], &info);
			printf("%-15lld\t%s\t%s\n", values[i],info.symbol,info.long_descr);
		}
		free(events);
		free(values);	
	}
}
