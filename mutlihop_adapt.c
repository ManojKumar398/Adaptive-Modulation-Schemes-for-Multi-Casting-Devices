#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include "nakagami_markov_chains.h"
#include "PER_values.h"

double run_session(int K, int info_bits_per_pkt, int num_first_hop_nodes, int num_second_hop_nodes, double MENR0, int m, int fdTs, int num_channel_states, int starting_idx, int adapt_type);
int configure_network(int num_first_hop_nodes, int num_second_hop_nodes, double MENR0, double nakagami_m, int num_channel_states, int starting_idx);
int markov_chain(int input_state, int m, int fdTs);
double compute_inst_MENR(int fading_state, int nominal_MENR, int m, int fdTs);
int packet_reception(double MENR, int MCS_idx);
int find_suggested_mcs(double MENR);
int adapt_mcs(int source_idx, int adapt_type, int current_mcs);
double compute_pkt_duration(int info_bits_per_pkt, int MCS_idx);
double bits_per_symbol[15] = {0.26,0.346,0.472,0.62,0.766,0.52,0.792,0.944,1.24,1.532,1.04,1.584,1.888,2.48,3.064};

__int64 x0longlong = 746384257, y0longlong = 344621271;
__int64 large_long;
double large_double, u0double, uni;
clock_t current;

typedef struct node_detail{
	int decoded; // 1 if node has decoded the file, 0 otherwise
	int pkts_received; // Number of packets received
	int pkts_sent;  // Number of packets sent
	double completion_time;
}node_details;

node_details *node;

typedef struct link_detail{
	int fading_state;
	double nominal_MENR;
	double MENR_offset;
	double instantaneous_MENR;
	int suggested_mcs;
	int next_hop;
}link_details;

link_details **link;

int num_nodes;

int main() {

	// Parameter values
	int K = 100; // Packets in the file
	int info_bits_per_pkt = 2400;
	int num_first_hop_nodes = 5; // No. of destinations
	int num_second_hop_nodes = 0; // Ignore for now
	double MENR0 = 5;
	int m = 250; // The Nakagami-m parameter multiplied by 100
	int fdTs = 20; // The normalized Doppler spread multiplied by 1000
	int num_channel_states = 12; // Number of states in the Markov chain
	int starting_idx = 14; // MCS index to be used in the beginning
	int adapt_type; // 0: No adaptation, 1: Min-rate, 2:Max-rate, 3: Mid-rate, 4: Max-DRR

	FILE *fp;

	int i;
	double session_throughput, avg_throughput;

	for (adapt_type = 0; adapt_type < 5; adapt_type++) {

		switch(adapt_type) {
			case 0:
				fp = fopen("data/fixed_rate.txt", "w");
				break;
			case 1:
				fp = fopen("data/min_rate.txt", "w");
				break;
			case 2:
				fp = fopen("data/max_rate.txt", "w");
				break;
			case 3:
				fp = fopen("data/mid_rate.txt", "w");
				break;
			case 4:
				fp = fopen("data/max_drr.txt", "w");
				break;
		}


		for (MENR0 = 0; MENR0 <= 20; MENR0 += 0.25) {
			avg_throughput = 0;
			for (i = 0; i < 1000; i++){
				session_throughput = run_session(K, info_bits_per_pkt, num_first_hop_nodes, num_second_hop_nodes, MENR0, m, fdTs, num_channel_states, starting_idx, adapt_type);
				//printf("session_throughput: %f\n", session_throughput);
				avg_throughput += session_throughput;
				//printf("avg_throughput: %f\n", avg_throughput/(i + 1));
			}
			printf("Nominal MENR: %f \t avg_throughput: %f\n", MENR0, avg_throughput/(i + 1));
			fprintf(fp, "%f\t %f\n", MENR0, avg_throughput/(i + 1));
		}

		fclose(fp);
	}


	/* fp = fopen(filename, "a");
	fprintf(fp, "%.2Lf\t%.6Lf\n", (double) snr, FER);
	fclose(fp); */

	return 0;
}


double run_session(int K, int info_bits_per_pkt, int num_first_hop_nodes, int num_second_hop_nodes, double MENR0, int m, int fdTs, int num_channel_states, int starting_idx, int adapt_type) {

	int MCS_idx, dest_idx, num_first_hop_decoded, num_tx = 0;
	double packet_duration, session_completion_time = 0, session_throughput;

	num_nodes = configure_network(num_first_hop_nodes, num_second_hop_nodes, MENR0, m, num_channel_states, starting_idx);
	//printf("Number of nodes (including the source): %d\n",num_nodes);

	MCS_idx = starting_idx;

	num_first_hop_decoded = 0;

	while (num_first_hop_decoded == 0 && num_tx < K*200) { // Continue as long as all first hop nodes are decoded (or too many transmissions are made and decoding still not complete)

		MCS_idx = adapt_mcs(0, adapt_type, MCS_idx);
		packet_duration = compute_pkt_duration(info_bits_per_pkt, MCS_idx);


		for (dest_idx = 1; dest_idx <= num_first_hop_nodes; dest_idx++)	{
			if (node[dest_idx].decoded == 0) {

				link[0][dest_idx].fading_state = markov_chain(link[0][dest_idx].fading_state, m, fdTs);
				link[0][dest_idx].instantaneous_MENR = compute_inst_MENR(link[0][dest_idx].fading_state, link[0][dest_idx].nominal_MENR, m, fdTs);
				link[0][dest_idx].suggested_mcs = find_suggested_mcs(link[0][dest_idx].instantaneous_MENR);
				//printf("MENR: %f \t MCS: %d\n",link[0][dest_idx].instantaneous_MENR, link[0][dest_idx].suggested_mcs);
				node[dest_idx].pkts_received +=  packet_reception(link[0][dest_idx].instantaneous_MENR, MCS_idx);

				if (node[dest_idx].pkts_received >= K) {
					node[dest_idx].decoded = 1;
					num_first_hop_decoded++;
				}

				node[dest_idx].completion_time = node[dest_idx].completion_time + packet_duration;
			}
		}

		num_tx++;
	}

	session_completion_time = 0;
	for (dest_idx = 1; dest_idx < num_nodes; dest_idx++)	{
		//printf("completion_time: %f \n", node[dest_idx].completion_time);
		if (node[dest_idx].completion_time > session_completion_time) {
			session_completion_time = node[dest_idx].completion_time;
		}
	}

	session_throughput = (double) K * info_bits_per_pkt /session_completion_time;

	//printf("session_throughput: %f \n", session_throughput);

	return session_throughput;
}



void uniformRV()
{
	large_long = 343597;
	large_long = large_long * 100000 + 38267;
	large_double = 171798;
	large_double = large_double * 100000 + 69184;
	x0longlong += time(NULL);
	current = clock();
	srand((unsigned int)(current*100));

	x0longlong=(x0longlong*8404997 + 1)&large_long;
	u0double=(long double)x0longlong/large_double;
	//printf("%d\n", u0double);
}


int configure_network(int num_first_hop_nodes, int num_second_hop_nodes, double MENR0, double m, int num_channel_states, int starting_idx) {

	int i, j, k;


	num_nodes = num_first_hop_nodes + num_second_hop_nodes + 1; // No. of nodes (including the source)

	//--------------------- Allocate Memory --------------------------//

	node = (node_details *) malloc(sizeof(node_details) * num_nodes);

	link = (link_details **) malloc(num_nodes*sizeof(link_details *));
	for(i = 0; i < num_nodes; i++) {link[i] = (link_details *) malloc(num_nodes * sizeof(link_details));}


	//--------------------- Initialize node status --------------------------//

	for (i = 0; i < num_nodes; i++) {
		node[i].pkts_received = 0;
		node[i].pkts_sent = 0;
		node[i].decoded = 0;
		node[i].completion_time = 0;
	}


	//-------------- Assign nominal link conditions-------------------//

	for (i = 0; i < num_nodes; i++) {
		for (j = 0; j < num_nodes; j++) {
			link[i][j].nominal_MENR = -100;  // Initialization
			link[i][j].MENR_offset = 0;
			link[i][j].next_hop = 0;
		}
	}

	link[0][1].MENR_offset = -10;
	link[0][2].MENR_offset = -10;
	link[0][3].MENR_offset = 0;
	link[0][4].MENR_offset = 0;
	link[0][5].MENR_offset = 0;

	// Links from the source to the first-hop nodes
	for (j = 1; j <= num_first_hop_nodes; j++) {
		link[0][j].nominal_MENR = MENR0  + link[0][j].MENR_offset;
		link[0][j].next_hop = 1;
	}

	// For now, suppose node 1 is the relay
	for (j = num_first_hop_nodes + 1; j < num_nodes; j++) {
		link[1][j].nominal_MENR = MENR0 + link[1][j].MENR_offset;
		link[1][j].next_hop = 1;
	}



	//-------------- Intialize link states and starting MCS -------------------//

	for (i = 0; i < num_nodes; i++) {
		for (j = 0; j < num_nodes; j++) {

			uniformRV();
			k = 0;
			while (num_channel_states * u0double > 2.0*k) {
				link[i][j].fading_state = k;
				k++;
			}

			link[i][j].suggested_mcs =  starting_idx;
		}
	}


	return num_nodes;
}



/*  Nakagami-m fading Markov Chain     */
int markov_chain(int input_state, int m, int fdTs)
{
	int output_state=0;
	long double u1;
	double temp;
	int i,count;

	uniformRV();u1 = u0double/2;

	switch (m) {
	case 100: // m = 1
		switch (fdTs) {
		case 200:  // fdTs = 0.200
			temp = p_matrix_m1_2[output_state][input_state];
			while (u1 > temp) {
				output_state++;
				temp = temp + p_matrix_m1_2[output_state][input_state];
			}
			break;
		case 100:  // fdTs = 0.100
			temp = p_matrix_m1_1[output_state][input_state];
			while (u1 > temp) {
				output_state++;
				temp = temp + p_matrix_m1_1[output_state][input_state];
			}
			break;
		case 50:  // fdTs = 0.050
			temp = p_matrix_m1_05[output_state][input_state];
			while (u1 > temp) {
				output_state++;
				temp = temp + p_matrix_m1_05[output_state][input_state];
			}
			break;
		case 20:  // fdTs = 0.020
			temp = p_matrix_m1_02[output_state][input_state];
			while (u1 > temp) {
				output_state++;
				temp = temp + p_matrix_m1_02[output_state][input_state];
			}
			break;
		case 10:  // fdTs = 0.010
			temp = p_matrix_m1_01[output_state][input_state];
			while (u1 > temp) {
				output_state++;
				temp = temp + p_matrix_m1_01[output_state][input_state];
			}
			break;
		case 5:  // fdTs = 0.005
			temp = p_matrix_m1_005[output_state][input_state];
			while (u1 > temp) {
				output_state++;
				temp = temp + p_matrix_m1_005[output_state][input_state];
			}
			break;
		case 2:  // fdTs = 0.002
			temp = p_matrix_m1_002[output_state][input_state];
			while (u1 > temp) {
				output_state++;
				temp = temp + p_matrix_m1_002[output_state][input_state];
			}
			break;
		}
		break;
	case 75: // m = 0.75
		switch (fdTs) {
		case 50:  // fdTs = 0.05
			temp = p_matrix_m075_05[output_state][input_state];
			while (u1 > temp) {
				output_state++;
				temp = temp + p_matrix_m075_05[output_state][input_state];
			}
			break;
		case 20:  // fdTs = 0.02
			temp = p_matrix_m075_02[output_state][input_state];
			while (u1 > temp) {
				output_state++;
				temp = temp + p_matrix_m075_02[output_state][input_state];
			}
			break;
		case 5:  // fdTs = 0.005
			temp = p_matrix_m075_005[output_state][input_state];
			while (u1 > temp) {
				output_state++;
				temp = temp + p_matrix_m075_005[output_state][input_state];
			}
			break;
		}
		break;
	case 90: // m = 0.9

		switch (fdTs) {
		case 50:  // fdTs = 0.05
			temp = p_matrix_m09_05[output_state][input_state];
			while (u1 > temp) {
				output_state++;
				temp = temp + p_matrix_m09_05[output_state][input_state];
			}
			break;
		case 20:  // fdTs = 0.02
			temp = p_matrix_m09_02[output_state][input_state];
			while (u1 > temp) {
				output_state++;
				temp = temp + p_matrix_m09_02[output_state][input_state];
			}
			break;
		case 5:  // fdTs = 0.005
			temp = p_matrix_m09_005[output_state][input_state];
			while (u1 > temp) {
				output_state++;
				temp = temp + p_matrix_m09_005[output_state][input_state];
			}
			break;
		}
		break;
	case 150: // m = 1.5
		switch (fdTs) {
		case 50:  // fdTs = 0.05
			temp = p_matrix_m15_05[output_state][input_state];
			while (u1 > temp) {
				output_state++;
				temp = temp + p_matrix_m15_05[output_state][input_state];
			}
			break;
		case 20:  // fdTs = 0.02
			temp = p_matrix_m15_02[output_state][input_state];
			while (u1 > temp) {
				output_state++;
				temp = temp + p_matrix_m15_02[output_state][input_state];
			}
			break;
		case 10:  // fdTs = 0.01
			temp = p_matrix_m15_01[output_state][input_state];
			while (u1 > temp) {
				output_state++;
				temp = temp + p_matrix_m15_01[output_state][input_state];
			}
			break;
		case 5:  // fdTs = 0.005
			temp = p_matrix_m15_005[output_state][input_state];
			while (u1 > temp) {
				output_state++;
				temp = temp + p_matrix_m15_005[output_state][input_state];
			}
			break;
		}
		break;
	case 250: // m = 2.5
		switch (fdTs) {
		case 50:  // fdTs = 0.05
			temp = p_matrix_m25_05[output_state][input_state];
			while (u1 > temp) {
				output_state++;
				temp = temp + p_matrix_m25_05[output_state][input_state];
			}
			break;
		case 20:  // fdTs = 0.02
			temp = p_matrix_m25_02[output_state][input_state];
			while (u1 > temp) {
				output_state++;
				temp = temp + p_matrix_m25_02[output_state][input_state];
			}
			break;
		case 10:  // fdTs = 0.01
			temp = p_matrix_m25_01[output_state][input_state];
			while (u1 > temp) {
				output_state++;
				temp = temp + p_matrix_m25_01[output_state][input_state];
			}
			break;
		case 5:  // fdTs = 0.005
			temp = p_matrix_m25_005[output_state][input_state];
			while (u1 > temp) {
				output_state++;
				temp = temp + p_matrix_m25_005[output_state][input_state];
			}
			break;
		}
		break;
	case 325: // m = 3.25
		switch (fdTs) {
		case 50:  // fdTs = 0.05
			temp = p_matrix_m325_05[output_state][input_state];
			while (u1 > temp) {
				output_state++;
				temp = temp + p_matrix_m325_05[output_state][input_state];
			}
			break;
		case 20:  // fdTs = 0.02
			temp = p_matrix_m325_02[output_state][input_state];
			while (u1 > temp) {
				output_state++;
				temp = temp + p_matrix_m325_02[output_state][input_state];
			}
			break;
		case 5:  // fdTs = 0.005
			temp = p_matrix_m325_005[output_state][input_state];
			while (u1 > temp) {
				output_state++;
				temp = temp + p_matrix_m325_005[output_state][input_state];
			}
			break;
		}
		break;
	case 576: // m = 5.76
		switch (fdTs) {
		case 50:  // fdTs = 0.05
			temp = p_matrix_m576_05[output_state][input_state];
			while (u1 > temp) {
				output_state++;
				temp = temp + p_matrix_m576_05[output_state][input_state];
			}
			break;
		case 20:  // fdTs = 0.02
			temp = p_matrix_m576_02[output_state][input_state];
			while (u1 > temp) {
				output_state++;
				temp = temp + p_matrix_m576_02[output_state][input_state];
			}
			break;
		case 5:  // fdTs = 0.005
			temp = p_matrix_m576_005[output_state][input_state];
			while (u1 > temp) {
				output_state++;
				temp = temp + p_matrix_m576_005[output_state][input_state];
			}
			break;
		}
		break;
	}

	return output_state;
}


int adapt_mcs(int source_idx, int adapt_type, int current_mcs) {

	int i, dest_idx, mcs_idx, max_drr_mcs, lowest_rate_mcs, highest_rate_mcs;
	double lowest_info_rate, highest_info_rate, mid_info_rate, data_recovery_rate, max_drr, diff;


	if (adapt_type == 0) {
		// No adaptation
		mcs_idx = current_mcs;
	}
	else if (adapt_type == 1) {
		// The min-rate protocol: Use the lowest-rate MCS among all the suggested MCSs
		mcs_idx = 14;
		for (dest_idx = 1; dest_idx < num_nodes; dest_idx++) {
			if (link[source_idx][dest_idx].next_hop &&  !node[dest_idx].decoded &&  bits_per_symbol[link[source_idx][dest_idx].suggested_mcs] < bits_per_symbol[mcs_idx]) {
				mcs_idx = link[source_idx][dest_idx].suggested_mcs;
			}
		}
	}
	else if (adapt_type == 2) {
		// The max-rate protocol: Use the highest-rate MCS among all the suggested MCSs
		mcs_idx = 0;
		for (dest_idx = 1; dest_idx < num_nodes; dest_idx++) {
			if (link[source_idx][dest_idx].next_hop &&  !node[dest_idx].decoded &&  bits_per_symbol[link[source_idx][dest_idx].suggested_mcs] > bits_per_symbol[mcs_idx]) {
				mcs_idx = link[source_idx][dest_idx].suggested_mcs;
			}
		}
	}
	else if (adapt_type == 3) {
		// The mid-rate protocol: Use the MCS whose rate is half-way between the lowest and the highest suggested MCSs
		lowest_rate_mcs = 14;
		highest_rate_mcs = 0;
		for (dest_idx = 1; dest_idx < num_nodes; dest_idx++) {
			if (link[source_idx][dest_idx].next_hop &&  !node[dest_idx].decoded &&  bits_per_symbol[link[source_idx][dest_idx].suggested_mcs] < bits_per_symbol[lowest_rate_mcs]) {
				lowest_rate_mcs = link[source_idx][dest_idx].suggested_mcs;
			}
			if (link[source_idx][dest_idx].next_hop &&  !node[dest_idx].decoded &&  bits_per_symbol[link[source_idx][dest_idx].suggested_mcs] > bits_per_symbol[highest_rate_mcs]) {
				highest_rate_mcs = link[source_idx][dest_idx].suggested_mcs;
			}
		}
		lowest_info_rate = bits_per_symbol[lowest_rate_mcs];
		highest_info_rate = bits_per_symbol[highest_rate_mcs];
		mid_info_rate = 0.5 * (lowest_info_rate + mid_info_rate);

		diff = 1000.0;
		for (i = 0; i < 15; i++) {
			if (fabs(bits_per_symbol[i] - mid_info_rate) < diff) {
				mcs_idx = i;
			}
		}
	}
	else if (adapt_type == 4) {
		// Max-DRR: Use the MCS that maximizes the cumulative data-recovery rate
		max_drr = 0;
		mcs_idx = 0;

		for (i = 0; i < 15; i++) {
			data_recovery_rate = 0;
			for (dest_idx = 1; dest_idx < num_nodes; dest_idx++) {
				if (link[source_idx][dest_idx].next_hop &&  !node[dest_idx].decoded) {
					data_recovery_rate += (bits_per_symbol[i] <= bits_per_symbol[link[source_idx][dest_idx].suggested_mcs]) * bits_per_symbol[i];
				}
			}
			if (data_recovery_rate > max_drr) {
				max_drr = data_recovery_rate;
				mcs_idx = i;
			}
		}
	}

	return mcs_idx;
}

int packet_reception(double MENR, int MCS_idx) {

	double error_prob;
	int MENR_idx;



	if (MENR < -6) {
		error_prob = 1;
	}
	else if (MENR > 13) {
		error_prob = 0;
	}
	else {
		MENR_idx = (int) ((MENR + 6.0)/0.25);
		error_prob = PER[MENR_idx][MCS_idx];
	}

	uniformRV();
	if (u0double/2.0 < error_prob) {
		return 0;
	}
	else {
		return 1;
	}
}


double compute_inst_MENR(int fading_state, int nominal_MENR, int m, int fdTs)
{
	double loss;
	double fading_delta = 2;

	switch (m) {

		case 100: // m = 1
			fading_delta = 2;
			loss = -fading_delta*(long double)fading_state + 16;
			break;
		case 75: // m = 0.75
			fading_delta = 2.5;
			loss = -fading_delta*(long double)fading_state + 20.75;
			break;
		case 90: // m = 0.9
			fading_delta = 2.25;
			loss = -fading_delta*(long double)fading_state + 18;
			break;
		case 150: // m = 1.5
			fading_delta = 1.75;
			loss = -fading_delta*(long double)fading_state + 13.5;
			break;
		case 250: // m = 2.5
			fading_delta = 1.25;
			loss = -fading_delta*(long double)fading_state + 9;
			break;
		case 325: // m = 3.25
			fading_delta = 1;
			loss = -fading_delta*(long double)fading_state + 7;
			break;
		case 576: // m = 5.76
			fading_delta = 0.75;
			loss = -fading_delta*(long double)fading_state + 5;
			break;
		default: // Static channel is default
			loss = 0;
			break;
	}

	loss = (fading_delta != 0) * loss;


	return nominal_MENR - loss;
}


double compute_pkt_duration(int info_bits_per_pkt, int MCS_idx) {

	double pkt_duration;

	//printf("info_bits_per_pkt: %d \t MCS_idx: %d\n",info_bits_per_pkt, MCS_idx);
	pkt_duration = info_bits_per_pkt/bits_per_symbol[MCS_idx];
	return pkt_duration;
}


int find_suggested_mcs(double MENR) {

	double per_symbol_throughput[15];
	int MCS_idx, MENR_idx, suggested_MCS_idx;


	if (MENR < -6) {
		for (MCS_idx = 0; MCS_idx < 15; MCS_idx++) {
			per_symbol_throughput[MCS_idx] = 0;
		}
	}
	else if (MENR > 13) {
		for (MCS_idx = 0; MCS_idx < 15; MCS_idx++) {
			per_symbol_throughput[MCS_idx] = bits_per_symbol[MCS_idx];
		}
	}
	else {
		MENR_idx = (int) ((MENR + 6.0)/0.25);
		for (MCS_idx = 0; MCS_idx < 15; MCS_idx++) {
			per_symbol_throughput[MCS_idx] = (1 - PER[MENR_idx][MCS_idx]) * bits_per_symbol[MCS_idx];
		}
	}

	//printf("MENR: %f\n", MENR);


	suggested_MCS_idx = 0;
	for (MCS_idx = 0; MCS_idx < 15; MCS_idx++) {
		//printf("mcs_idx: %d \t per_symbol_throughput: %f", MCS_idx, per_symbol_throughput[MCS_idx]);
		if (per_symbol_throughput[MCS_idx]  >= per_symbol_throughput[suggested_MCS_idx]) {
			suggested_MCS_idx = MCS_idx;
		}
	}

	return suggested_MCS_idx;
}
