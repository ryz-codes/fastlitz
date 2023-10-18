#ifndef KPU_2_2_H
#define KPU_2_2_H

#define Q_NUM_WEIGHTS 	5 
#define Q_NUM_DIMS 	2 

static const double Qweights[5] = 
	{0.2777778,
	 0.2777778,
	 -0.1111111999999999,
	 0.2777778,
	 0.2777778};

static const double Qnodes[2][5] = {
	{0.1127017,
	 0.5,
	 0.5,
	 0.5,
	 0.8872983},

	{0.5,
	 0.1127017,
	 0.5,
	 0.8872983,
	 0.5}
};

#endif /** define KPU_2_2_H **/
