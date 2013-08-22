#include "XSbench_header.h"

void counter_init( int *eventset, int *num_papi_events )
{
	char error_str[PAPI_MAX_STR_LEN];
	//  int events[] = {PAPI_TOT_INS,PAPI_BR_INS,PAPI_SR_INS};
	int events[] = {PAPI_TOT_INS,PAPI_LD_INS,PAPI_FP_INS};
	int stat;

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
