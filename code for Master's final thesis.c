# include "udf.h"
# include "stdio.h"

# define Oil_Domain 2
# define Gas_Domain 3
# define PVT_tables "C:/thesis/case/CFDdata/pvtvalsSPE1.txt"
# define initial_values "C:/thesis/case/CFDdata/inits.txt"
# define testouts "C:/thesis/case/CFDdata/busa.txt"
# define itersperts 20

#define cg 5.28e-5			/* gas compressibility (/psi) */
#define co 7.00e-7			/* oil compressibility (/psi) */
#define rho_w_init 1088.9			/* initial water density (kg/m3)*/
#define rho_o_init 855			/* initial oil density (kg/m3)*/
#define p_init 14.6959488			/* initial pressure (psi) */

#define psi_to_pascal 6894.757 		/* conversion factor from psi to Pascal */

#define porosity 0.2
#define equil_time 600


double curts;
double tfin;
double curdelta;

double oil_pvt[500][4];
double gas_pvt[50][3];
double oil_surf_density;
double gas_surf_density;
double init_Rs;
double init_gvf;
int oil_line_counter=0;
int gas_line_counter=0;

int mf = 0;

double init_gmf_in = 0.1;

double Pmin;
double Pmax;

void fillPminmax()
{
    int i = 0;
	double PmaxPrev = 0;
	Pmin = 0;
	Pmax = 0;
	for (i = 0; i < oil_line_counter; i++)
	{
		if (oil_pvt[i][0] != -1)
		{
			if (Pmin == 0)
			{
				Pmin = oil_pvt[i][1];
			}
			else if (Pmin > oil_pvt[i][1]) Pmin = oil_pvt[i][1];
			if (Pmax == 0)
			{
				Pmax = oil_pvt[i][1];
			}
			else if (Pmax < oil_pvt[i][1])
			{
				PmaxPrev = Pmax;
				Pmax = oil_pvt[i][1];
			}
		}
	}
	Pmax = PmaxPrev;
}




typedef struct rsouts {
	double Bo;
	double ovisc;
} trsoutdat;


typedef struct pouts {
	double Rs;
	double P;
	double Bo;
	double ovisc;
	double Bg;
	double gvisc;
} tpoutdat;


typedef struct outs {
	double gvf;
	double gio;
	double od;
	double ov;
	double gd;
	double gv;
	double gsound;
	double osound;
	double mkrog;
	double mkrgo;
}toutdat;


trsoutdat singleRsInterpolate(double Rs)
{
	trsoutdat odat;


	int prev_rs_line = 0;
	int next_rs_line = 0;
	double prev_line[4];
	double next_line[4];
	double int_line[4];


	int f = 0;
	int i = 0;
	int j = 0;

	while (f == 0)
	{
		if (oil_pvt[i][0] != -1)
		{
			if (oil_pvt[i][0] <= Rs)
			{
				prev_rs_line = i;
				for (j = 0; j < 4; j++) prev_line[j] = oil_pvt[i][j];
				if (oil_pvt[i][0] == Rs) f = 1;
			}
			else if ((oil_pvt[i][0]>Rs) && (f != 1))
			{
				next_rs_line = i;
				for (j = 0; j < 4; j++) next_line[j] = oil_pvt[i][j];
				f = 1;
			}
		}
		i++;
	}

	if (next_rs_line != 0)
	{
		double prevRs = prev_line[0];
		double nextRs = next_line[0];
		int_line[0] = Rs;
		int_line[1] = (next_line[1] - prev_line[1]) / (nextRs - prevRs)*(Rs - prevRs) + prev_line[1];
		int_line[2] = (next_line[2] - prev_line[2]) / (nextRs - prevRs)*(Rs - prevRs) + prev_line[2];
		int_line[3] = (next_line[3] - prev_line[3]) / (nextRs - prevRs)*(Rs - prevRs) + prev_line[3];
	}
	else
	{
		for (i = 0; i < 4; i++) int_line[i] = prev_line[i];
	}

	odat.Bo = int_line[2];
	odat.ovisc = int_line[3];

	return(odat);
}

trsoutdat RsInterpolate(double Rs, double P)
{

	trsoutdat odat;

	int prev_rs_line = 0;
	int next_rs_line = 0;
	int prev_block_lc = 0;
	int next_block_lc = 0;
	double prev_block[50][4];
	double next_block[50][4];
	double int_block[50][4];
		int lc = 0;
	double BoPrev;
	double BoNext;
	double PPrev;
	double PNext;
	double viscPrev;
	double viscNext;
	int oil_interp_lc;

	int f = 0;
	int i = 0;
	int j = 0;
	int k = 0;

	while (f == 0)
	{
		if (oil_pvt[i][0] != -1)
		{
			if (oil_pvt[i][0] <= Rs)
			{
				prev_rs_line = i;
				prev_block_lc = 1;
				for (j = 0; j < 4; j++) prev_block[0][j] = oil_pvt[i][j];
				if (oil_pvt[i][0] == Rs) f = 1;
				i++;
				while ((i<oil_line_counter) && (oil_pvt[i][0] == -1))
				{
					for (j = 0; j<4; j++) prev_block[prev_block_lc][j] = oil_pvt[i][j];
					i++;
					prev_block_lc++;
				}
			}
			else if ((oil_pvt[i][0]>Rs) && (f != 1))
			{
				next_rs_line = i;
				next_block_lc = 1;
				for (j = 0; j < 4; j++) next_block[0][j] = oil_pvt[i][j];
				i++;
				while ((i<oil_line_counter) && (oil_pvt[i][0] == -1))
				{
					for (j = 0; j<4; j++) next_block[next_block_lc][j] = oil_pvt[i][j];
					i++;
					next_block_lc++;
				}
				f = 1;
			}
		}
		oil_interp_lc = prev_block_lc;
		if (next_block_lc == 0)
		{
			for (j = 0; j < oil_interp_lc; j++)
			{
				for (k = 0; k < 4; k++)
				{
					int_block[j][k] = prev_block[j][k];
				}
			}
		}
		else
		{
			double prevRs = prev_block[0][0];
			double nextRs = next_block[0][0];
			int_block[0][0] = Rs;
			int_block[0][1] = (next_block[0][1] - prev_block[0][1]) / (nextRs - prevRs)*(Rs - prevRs) + prev_block[0][1];
			int_block[0][2] = (next_block[0][2] - prev_block[0][2]) / (nextRs - prevRs)*(Rs - prevRs) + prev_block[0][2];
			int_block[0][3] = (next_block[0][3] - prev_block[0][3]) / (nextRs - prevRs)*(Rs - prevRs) + prev_block[0][3];
			for (j = 1; j < prev_block_lc; j++)
			{
				if (j != 0)
				{
					int_block[j][0] = prev_block[j][0];
					int_block[j][1] = prev_block[j][1];
				}
				int_block[j][2] = (next_block[j - 1][2] - prev_block[j][2]) / (nextRs - prevRs)*(Rs - prevRs) + prev_block[j][2];
				int_block[j][3] = (next_block[j - 1][3] - prev_block[j][3]) / (nextRs - prevRs)*(Rs - prevRs) + prev_block[j][3];
			}
		}
	}




		f = 0;
    if (int_block[lc][1] > P)
    {
        real deltaBo;
        real deltaVisc;        
        
        deltaBo=(int_block[1][2]-int_block[0][2])/(int_block[1][1]-int_block[0][1])*(int_block[0][1]-P);
        deltaVisc=(int_block[1][3]-int_block[0][3])/(int_block[1][1]-int_block[0][1])*(int_block[0][1]-P);
        
        odat.Bo = int_block[0][2]-deltaBo;
	    odat.ovisc = int_block[0][3]-deltaVisc;

       /* BoPrev = int_block[lc][2];
	    viscPrev = int_block[lc][3];
	    PPrev = 1;
        BoNext = int_block[lc][2];
	    viscNext = int_block[lc][3];
	    PNext = 1;*/
    }    
    else
    {
        while (f == 0)
	    {
		    if (int_block[lc][1] <= P)
		    {
			    BoPrev = int_block[lc][2];
			    viscPrev = int_block[lc][3];
			    PPrev = int_block[lc][1];
		    }
		    if (int_block[lc][1] > P)
		    {
			    BoNext = int_block[lc][2];
			    viscNext = int_block[lc][3];
			    PNext = int_block[lc][1];
			    f = 1;
		    }
		    lc++;
	    }
    	odat.Bo = (BoNext - BoPrev) / (PNext - PPrev)*(P - PPrev) + BoPrev;
	    odat.ovisc = (viscNext - viscPrev) / (PNext - PPrev)*(P - PPrev) + viscPrev;
    }

	return(odat);

}



tpoutdat Pinterpolate(double pressure)
{

	int gc = 0;
	int oc = 0;
	int f = 0;

	FILE *fil;

	tpoutdat odat;

	double oil_prev_line[4];
	double oil_next_line[4];
	double gas_prev_line[3];
	double gas_next_line[3];

	double P = pressure;

	int i = 0;

	if (P > Pmax) P = Pmax;

	while (f == 0)
	{
		if (oil_pvt[oc][0] != -1)
		{
			if (P > oil_pvt[oc][1])
			{
				for (i = 0; i < 4; i++) oil_prev_line[i] = oil_pvt[oc][i];
				for (i = 0; i < 3; i++) gas_prev_line[i] = gas_pvt[gc][i];
			}
			if (P <= oil_pvt[oc][1])
			{
				for (i = 0; i < 4; i++) oil_next_line[i] = oil_pvt[oc][i];
				for (i = 0; i < 3; i++) gas_next_line[i] = gas_pvt[gc][i];
				f = 1;
			}
			gc++;
		}
		oc++;
	}

	odat.P = pressure;
	odat.Rs = (oil_next_line[0] - oil_prev_line[0]) / (oil_next_line[1] - oil_prev_line[1])*(P - oil_prev_line[1]) + oil_prev_line[0];
	odat.Bo = (oil_next_line[2] - oil_prev_line[2]) / (oil_next_line[1] - oil_prev_line[1])*(P - oil_prev_line[1]) + oil_prev_line[2];
	odat.ovisc = (oil_next_line[3] - oil_prev_line[3]) / (oil_next_line[1] - oil_prev_line[1])*(P - oil_prev_line[1]) + oil_prev_line[3];
	odat.Bg = (gas_next_line[1] - gas_prev_line[1]) / (gas_next_line[0] - gas_prev_line[0])*(pressure - gas_prev_line[0]) + gas_prev_line[1];
	odat.gvisc = (gas_next_line[2] - gas_prev_line[2]) / (gas_next_line[0] - gas_prev_line[0])*(pressure - gas_prev_line[0]) + gas_prev_line[2];

	return(odat);
}



toutdat taking_pars(double pressure, double g_vf)
{

	toutdat dout;

	FILE *fil;

	double n_dens_oil;
	double n_dens_gas;
	double n_visc_oil;
	double n_visc_gas;
	double n_g_vf = g_vf;

	double Rs = init_Rs;

	double n_GiOmf = (Rs*gas_surf_density)/(oil_surf_density+Rs*gas_surf_density);


	int gc = 0;
	int oc = 0;
	int f = 0;

	tpoutdat pintout;
	trsoutdat rsintout;

	trsoutdat dprsout;
	tpoutdat dppout;
	double dp = 10;
	double dpodens = 0;
	double dpgdens = 0;



	pintout = Pinterpolate(pressure);
	rsintout = RsInterpolate(Rs, pressure);

	dppout = Pinterpolate(pressure + dp);
	dprsout = RsInterpolate(Rs, pressure + dp);

	dpodens = (oil_surf_density + gas_surf_density*Rs) / dprsout.Bo;
	dpgdens = gas_surf_density / dppout.Bg;

	n_dens_oil = (oil_surf_density + gas_surf_density*Rs) / rsintout.Bo;
	n_visc_oil = rsintout.ovisc;
	n_dens_gas = gas_surf_density / pintout.Bg;
	n_visc_gas = pintout.gvisc;

	dout.gsound = sqrt(1e5*dp / (dpgdens - n_dens_gas));
	dout.osound = sqrt(1e5*dp / (dpodens - n_dens_oil));
	
	dout.gvf = n_g_vf;
	dout.gio = n_GiOmf;
	dout.od = n_dens_oil;
	dout.ov = n_visc_oil;
	dout.gd = n_dens_gas;
	dout.gv = n_visc_gas;

	return(dout);
}


toutdat changing_pars(double pressure, double g_vf, double GiOmf, double Vcell, double dens_oil, double dens_gas, double visc_oil, double visc_gas)
{

	toutdat dout;

	toutdat pdout;

	tpoutdat pintout;
	trsoutdat rsintout;

	trsoutdat dprsout;
	tpoutdat dppout;
	trsoutdat sodat;
	double dp = 10;
	double dpodens = 0;
	double dpgdens = 0;
		double n_GiOmf = GiOmf;
		double n_dens_oil = dens_oil;
	double n_dens_gas = dens_gas;
	double n_visc_oil = visc_oil;
	double n_visc_gas = visc_gas;
	double n_g_vf = g_vf;

	int gc = 0;
	int oc = 0;
	int f = 0;


	int oil_interp_lc = 0;
	double oil_interp[50][4];

	double prev_block[50][4];
	double next_block[50][4];


	double ovf = 1.0 - g_vf;
	double n_ovf = ovf;
	double mtk;
	double R;
	double Rs_in;
	FILE *fil;

	pintout = Pinterpolate(pressure);
	dppout = Pinterpolate(pressure + dp);

	Rs_in = (oil_surf_density*GiOmf) / (gas_surf_density*(1 - GiOmf));

	sodat = singleRsInterpolate(Rs_in);
	
	/*dens_oil= (oil_surf_density + gas_surf_density*Rs_in) /sodat.Bo;
	dens_gas = gas_surf_density / pintout.Bg;*/


	
	if (ovf == 0.0)
	{
		n_GiOmf = 0;
		n_dens_gas = gas_surf_density / pintout.Bg;
		n_visc_gas = pintout.gvisc;
		n_visc_oil = pintout.ovisc;
		n_dens_oil = (oil_surf_density + gas_surf_density*pintout.Rs) / pintout.Bo;
		R = 0.0;
	}
	else
	{
		R = ((g_vf*dens_gas + ovf*dens_oil*GiOmf)*oil_surf_density) / (ovf*dens_oil*(1 - GiOmf)*gas_surf_density);
	}

	if ((R > pintout.Rs)&&(R!=0.0))
	{
		n_ovf = (ovf*dens_oil*(1 - GiOmf)*pintout.Bo) / oil_surf_density;
		if (n_ovf > 1.0) n_ovf = 1.0;
		n_dens_oil = (oil_surf_density + gas_surf_density*pintout.Rs) / pintout.Bo;
		n_GiOmf = (ovf*dens_oil*(1 - GiOmf)*pintout.Rs*gas_surf_density) / (n_ovf*n_dens_oil*oil_surf_density);
		n_visc_oil = pintout.ovisc;
		n_dens_gas = gas_surf_density / pintout.Bg;
		n_visc_gas = pintout.gvisc;
		dprsout = RsInterpolate(pintout.Rs, pressure + dp);
		dpodens = (oil_surf_density + gas_surf_density*pintout.Rs) / dprsout.Bo;
		dpgdens = gas_surf_density / dppout.Bg;
	}
	else if ((R < pintout.Rs)&&(R!=0.0))
	{	
	mtk= dens_oil*ovf - n_dens_oil*n_ovf;
		n_ovf = 1.0;
		n_GiOmf = (g_vf*dens_gas + ovf*dens_oil*GiOmf) / (g_vf*dens_gas + ovf*dens_oil);
		n_dens_gas = gas_surf_density / pintout.Bg;
		n_visc_gas = pintout.gvisc;

		rsintout = RsInterpolate(R,pressure);

		dprsout = RsInterpolate(R, pressure + dp);
		dpodens = (oil_surf_density + gas_surf_density*R) / dprsout.Bo;
		dpgdens = gas_surf_density / dppout.Bg;

		n_dens_oil = (oil_surf_density + gas_surf_density*R) / rsintout.Bo;
		n_visc_oil = rsintout.ovisc;
	}


	dout.gsound = sqrt(1e5*dp / (dpgdens - n_dens_gas));
	dout.osound = sqrt(1e5*dp / (dpodens - n_dens_oil));

	/*if ((1.0 - n_ovf) != 0.0) n_dens_gas = (dens_gas*(1 - ovf) + dens_oil*ovf - n_dens_oil*n_ovf) / (1.0 - n_ovf);*/


	dout.mkrgo = 0.0;
	dout.mkrog = 0.0;



	if (mtk > 0.0)
	{
		dout.mkrog = mtk/(ovf*CURRENT_TIMESTEP);
	}
	else if (mtk < 0.0)
	{
		dout.mkrgo = (dens_gas*(1 - ovf)-n_dens_gas*(1 - n_ovf))/ ((1-ovf)*curdelta);
	}

	dout.gvf = 1.0 - n_ovf;
	dout.gio = n_GiOmf;
	dout.od = n_dens_oil;
	dout.ov = n_visc_oil;
	dout.gd = n_dens_gas;
	dout.gv = n_visc_gas;

	/*if (n_ovf > ovf)
	{
		fil = fopen("c:\\CFDdata\\busa.txt", "a");
		fprintf(fil, "%7.4f : %7.4f : %7.4f \n", pressure, GiOmf, dens_oil);
		fclose(fil);
	}*/
	return(dout);
}

toutdat propcalc(double pressure, double GiOmf)
{	double dp = 10;
	double dpodens = 0;
	double dpgdens = 0;
double Rs_in;
	double n_dens_oil;
	double n_dens_gas;
	double n_visc_oil;
	double n_visc_gas;

	int gc = 0;
	int oc = 0;
	int f = 0;


	int oil_interp_lc = 0;
	double oil_interp[50][4];

	double prev_block[50][4];
	double next_block[50][4];
	toutdat dout;

	toutdat pdout;
	trsoutdat sodat;
	tpoutdat pintout;
	trsoutdat rsintout;

	trsoutdat dprsout;
	tpoutdat dppout;

	FILE *fil;

	pintout = Pinterpolate(pressure);
	dppout = Pinterpolate(pressure + dp);


	



	Rs_in = (oil_surf_density*GiOmf) / (gas_surf_density*(1 - GiOmf));


	/*dens_oil= (oil_surf_density + gas_surf_density*Rs_in) /sodat.Bo;
	dens_gas = gas_surf_density / pintout.Bg;*/




    n_dens_gas = gas_surf_density / pintout.Bg;
	n_visc_gas = pintout.gvisc;

	rsintout = RsInterpolate(Rs_in, pressure);

	dprsout = RsInterpolate(Rs_in, pressure + dp);
	dpodens = (oil_surf_density + gas_surf_density*Rs_in) / dprsout.Bo;
	dpgdens = gas_surf_density / dppout.Bg;

	n_dens_oil = (oil_surf_density + gas_surf_density*Rs_in) / rsintout.Bo;
	n_visc_oil = rsintout.ovisc;

	dout.gsound = sqrt(1e5*dp / (dpgdens - n_dens_gas));
	dout.osound = sqrt(1e5*dp / (dpodens - n_dens_oil));


	dout.mkrgo = 0.0;
	dout.mkrog = 0.0;

	/*double mtk = dens_oil*ovf - n_dens_oil*n_ovf;*/

	
	dout.od = n_dens_oil;
	dout.ov = n_visc_oil;
	dout.gd = n_dens_gas;
	dout.gv = n_visc_gas;

	/*if (n_ovf > ovf)
	{
	fil = fopen("c:\\CFDdata\\busa.txt", "a");
	fprintf(fil, "%7.4f : %7.4f : %7.4f \n", pressure, GiOmf, dens_oil);
	fclose(fil);
	}*/
	return(dout);
}

DEFINE_PROPERTY(oil_speed_of_sound, cell, thread)
{	toutdat dout;
    real ao;
	real gio;
	real pressure;

	gio = C_UDSI(cell, thread, 0);
	pressure = C_P(cell, thread)+101325;



	dout = propcalc(pressure / 100000, gio);

	ao = dout.osound;

	/*ao = C_UDMI(cell, thread, 10);*/

	return(ao);
}

DEFINE_PROPERTY(gas_speed_of_sound, cell, thread)
{
    real ag;



	ag = 1000;

	/*ag = C_UDMI(cell, thread, 11);*/

	return(ag);
}


DEFINE_PROPERTY(set_oil_density, cell, thread)
{
    real ro;

	real gio;
	real pressure;

	FILE *fil;

	toutdat dout;

	gio = C_UDSI(cell, thread, 0);
	pressure = C_P(cell, thread)+101325;

	dout = propcalc(pressure / 100000, gio);

	ro = dout.od;
	
	/*ro = C_UDMI(cell, thread, 4);*/

	return(ro);
}

DEFINE_PROPERTY(set_gas_density, cell, thread)
{
    real rg;


	rg=123;

	/*rg = C_UDMI(cell, thread, 6);*/

	return(rg);
}

DEFINE_PROPERTY(set_oil_viscosity, cell, thread)
{
    real mo;

	real gio;
	real pressure;

	toutdat dout;


	gio = C_UDSI(cell, thread, 0);
	pressure = C_P(cell, thread)+101325;

	dout = propcalc(pressure / 100000, gio);

	mo = dout.ov*1e-3;

	/*mo = C_UDMI(cell, thread, 12)*1e-3;*/

	return(mo);
}

DEFINE_PROPERTY(set_gas_viscosity, cell, thread)
{
    real mg;

	real gio;
	real pressure;

	toutdat dout;

	Domain *odomain;
	Thread *othread;
	cell_t ocell;

	odomain = Get_Domain(Oil_Domain);

	thread_loop_c(othread, odomain)
	{
		begin_c_loop(ocell, othread)
		{
			if (ocell == cell)
			{
				gio = C_UDSI(ocell, othread, 0);
				pressure = C_P(ocell, othread)+101325;
			}
		}
		end_c_loop(ocell, othread)
	}

	dout = propcalc(pressure/100000, gio);

	mg = dout.gv*1e-3;

	/*mg = C_UDMI(cell, thread, 13)*1e-3;*/

	return(mg);
}

DEFINE_INIT(init_sequence, domain)
{
	Domain *odomain;
	Domain *gdomain;
	cell_t c;
	Thread *t;

	toutdat odat;

	FILE *fil;

    /*Reading the initial values*/

	fil = fopen(initial_values, "r");
	fscanf(fil, "%le", &init_Rs);
	fscanf(fil, "%le", &init_gvf);
	fclose(fil);
	init_Rs=332.3487;
	init_gvf=0.0;
	odomain = Get_Domain(Oil_Domain);
	gdomain = Get_Domain(Gas_Domain);

    curts=0.0;
    curdelta=0.0;
    tfin=0.0;

	thread_loop_c(t, odomain)
	{
		begin_c_loop(c, t)
		{
			odat = taking_pars((C_P(c, t)+101325) / 100000, 1.0 - init_gvf);
			C_VOF(c, t) = 1.0 - init_gvf;
			C_UDMI(c, t, 0) = C_P(c, t)+101325;
			C_UDMI(c, t, 1) = odat.gio;
			C_UDMI(c, t, 2) = odat.gio;
            C_UDMI(c,t,3)=odat.od;
            C_UDMI(c,t,4)=odat.od;
            C_UDMI(c,t,5)=odat.gd;
            C_UDMI(c,t,6)=odat.gd;
            C_UDMI(c,t,7)=init_gvf;
            C_UDMI(c,t,8)=init_gvf;
			C_UDMI(c, t, 10) = odat.osound;
			C_UDMI(c, t, 11) = odat.gsound;
			C_UDMI(c, t, 12) = odat.ov;
			C_UDMI(c, t, 13) = odat.gv;

			C_UDMI(c, t, 14) = 0.0;
			C_UDMI(c, t, 15) = 0.0;
			C_UDMI(c, t, 16) = 0.0;
			C_UDMI(c, t, 17) = 0.0;
			C_UDMI(c, t, 18) = 0.0;


			C_R(c, t) = odat.od;
			C_MU_L(c, t) = odat.ov*1e-3;
			C_UDSI(c, t, 0) = odat.gio;
			init_gmf_in = odat.gio;
		}
		end_c_loop(c, t)
	}
	thread_loop_c(t, gdomain)
	{
		begin_c_loop(c, t)
		{
			C_VOF(c, t) = init_gvf;
			C_R(c, t) = C_UDMI(c, t, 5);
			C_MU_L(c, t) = C_UDMI(c, t, 9)*1e-3;
		}
		end_c_loop(c, t)
	}
	Message("Solution is initialized \n");
	Message("Equilibration time - ", equil_time, "\n");
}


DEFINE_EXECUTE_AT_END(PhaseCorrection)
{
	Domain *domain;
	Domain *gdomain;
	Domain *mixdomain;
	cell_t c;
	cell_t gc;
	face_t f;
	Thread *t;
	Thread *gpt;
	int n;
 	int f_n;   
	face_t cf;
	Thread *ft;

	FILE *fil;

	toutdat datout;

	double press_in;
	double gvf_in;
	double gio_in;
	double cv_in;
	double odens_in;
	double ovisc_in;
	double gdens_in;
	double gvisc_in;

	int i = 0;  
	curts=CURRENT_TIMESTEP;
    tfin=curts+CURRENT_TIME;
    curdelta=curts;

	domain = Get_Domain(Oil_Domain);
	gdomain = Get_Domain(Gas_Domain);
	mixdomain = Get_Domain(1);




/*
	thread_loop_c(t, domain)
	{
		begin_c_loop(c, t)
		{
			C_R(c, t) = C_UDMI(c, t, 4);

			thread_loop_f(ft, domain)
			{
				begin_f_loop(cf, ft)
				{
					if (F_C0(cf, ft) == c)
					{

						cell_t adc;

						real flf;
						real dgio;

						flf = F_FLUX(cf, ft);

						if ((flf < 0) && (!BOUNDARY_FACE_THREAD_P(ft)))
						{
							adc = F_C1(cf, ft);
							dgio = C_UDSI(adc, t, 0);
						}
						else if ((flf < 0) && (BOUNDARY_FACE_THREAD_P(ft)))
						{
							dgio = init_gmf_in;
						}

						C_UDMI(c, t, 1) = dgio;

					}
				}
				end_f_loop(cf, ft)
			}

		}
		end_c_loop(c, t)
	}

*/

	thread_loop_c(t, mixdomain)
	{
		begin_c_loop(c, t)
		{
			press_in = C_P(c, t)+101325;
			C_UDMI(c, t, 0) = press_in / 100000;
		}
		end_c_loop(c, t)
	}

	thread_loop_c(t, domain)
	{
		begin_c_loop(c, t)
		{
			C_UDMI(c, t, 3) = C_R(c, t);
			C_UDMI(c, t, 9) = C_MU_L(c, t);
 		}
		end_c_loop(c, t)
	}

  thread_loop_c (gpt,gdomain)
    {
		begin_c_loop(gc, gpt)
		{
			/*C_UDMI(gc, gpt, 1) = C_UDMI(gc, gpt, 2);*/
            press_in = C_UDMI(gc, gpt, 0);
			gvf_in = C_VOF(gc, gpt);
			gio_in = C_UDMI(gc, gpt, 1);
			cv_in = C_VOLUME(gc, gpt);
			odens_in = C_UDMI(gc, gpt, 3);
			ovisc_in = C_UDMI(gc, gpt, 9);
			gdens_in = C_R(gc, gpt);
			gvisc_in = C_MU_L(gc, gpt);
			C_UDMI(gc, gpt, 5) = gdens_in;
			C_UDMI(gc, gpt, 7) = gvf_in;
			datout = changing_pars(press_in, gvf_in, gio_in, cv_in, odens_in, gdens_in, ovisc_in, gvisc_in);
			C_UDMI(gc, gpt, 8) = datout.gvf;
			C_UDMI(gc, gpt, 2) = datout.gio;
			C_UDMI(gc, gpt, 4) = datout.od;
//			C_UDMI(gc, gpt, 4) = 855;
			C_UDMI(gc, gpt, 6) = datout.gd;
//			C_UDMI(gc, gpt, 6) = 1.293;
			C_R(gc, gpt) = datout.gd;
//			C_R(gc, gpt) = 1.293;
			C_UDMI(gc, gpt, 10) = datout.osound;
			C_UDMI(gc, gpt, 11) = datout.gsound;
			C_UDMI(gc, gpt, 12) = datout.ov;
			C_UDMI(gc, gpt, 13) = datout.gv;
			C_UDMI(gc, gpt, 14) = odens_in*C_VOLUME(gc, gpt)*porosity*(1 - gvf_in)*gio_in;
			C_UDMI(gc, gpt, 15) = datout.od*C_VOLUME(gc, gpt)*porosity*(1 - datout.gvf)*datout.gio;
			C_UDMI(gc, gpt, 16) = odens_in*C_VOLUME(gc, gpt)*porosity*(1 - gvf_in);
						
		}
       end_c_loop (gc,gpt)
    }
  
  thread_loop_c(t, domain)
  {
	  begin_c_loop(c, t)
	  {
		  


		  real totfluxg = 0.0;
		  real totfluxo = 0.0;

		  real mtc = 0.0;
		  real giotc = 0.0;

		  real dgio = C_UDMI(c, t, 2);
C_R(c, t) = C_UDMI(c, t, 4);
		  thread_loop_f(ft, domain)
		  {
			  begin_f_loop(cf, ft)
			  {
				  if (F_C0(cf, ft) == c)
				  {
					  real flf;

					  cell_t adc;

					  flf = F_FLUX(cf, ft)*CURRENT_TIMESTEP;

					  if ((flf < 0) && (!BOUNDARY_FACE_THREAD_P(ft)))
					  {
						  adc = F_C1(cf, ft);
						  dgio = C_UDSI(adc, t, 0);
						  /*F_UDSI(cf, ft, 0) = dgio;*/
					  }
					  else if ((flf < 0) && (BOUNDARY_FACE_THREAD_P(ft)))
					  {
						  dgio = init_gmf_in;
					  }


//					  totfluxg += dgio*flf;
					  totfluxo += flf;
				  }
			  }
			  end_f_loop(cf, ft)
		  }

		  if (C_UDMI(c, t, 7) == C_UDMI(c, t, 8))
		  {
			  mtc = 0.0;
			  giotc = 0.0;
		  }
		  else
		  {

			 

			  real nmo;
			  real nmg;
			  real ngio;

 mtc = (C_UDMI(c, t, 15) - C_UDMI(c, t, 14) + totfluxg) / (C_VOLUME(c, t)*porosity);
			  nmo = C_UDMI(c, t, 16) - totfluxo;
			  nmg = C_UDMI(c, t, 14) - totfluxg;

			  ngio = nmg / nmo;


			  giotc = (C_UDMI(c, t, 2) - ngio);
		  }
		  C_UDMI(c, t, 17) = mtc;
		  C_UDMI(c, t, 18) = giotc;



	  }
	  end_c_loop(c, t)
  }

}


DEFINE_PROFILE(gmf_inlet_profile, t, i)
{
	face_t f;
	begin_f_loop(f, t)
	{
		F_PROFILE(f, t, i) = init_gmf_in;
	}
	end_f_loop(f,t)
}

DEFINE_PROFILE(uds_profile, t, i)
{
	cell_t c;
	begin_c_loop(c, t)
	{
        real uds1;
        real uds2;
        real uds_cur;
        real num_iter;		
        real itnum=itersperts;
        if (itnum!=0)
        {        
            uds1=C_UDMI(c,t,1);
            uds2=C_UDMI(c,t,2);
            num_iter=C_UDMI(c,t,17);
            uds_cur=uds1+num_iter*(uds2-uds1)/itnum;
        }
        else uds_cur=C_UDMI(c,t,1);
        C_PROFILE(c, t, i) = uds_cur;
	}
	end_c_loop(c, t)
}

DEFINE_PROFILE(gmf_outlet_profile, t, i)
{
	face_t f;
	cell_t c;
	Thread *ct;

	begin_f_loop(f, t)
	{		double gvf;
		c = F_C0(f, t);
		ct = THREAD_T0(t);

		gvf = C_UDSI(c, ct, 0);
		F_PROFILE(f, t, i) = gvf;
	}
	end_f_loop(f, t)
}

DEFINE_ADJUST(adjustment, dom)
{
	Domain *domain;
	Domain *gdomain;
	Domain *mixdomain;
	cell_t c;
	cell_t gc;
	face_t f;
	Thread *t;
	Thread *gpt;
	int n;
	FILE *fil;

	toutdat datout;

	double press_in;
	double gvf_in;
	double gio_in;
	double cv_in;
	double odens_in;
	double ovisc_in;
	double gdens_in;
	double gvisc_in;

	int i = 0;
	if (curts != CURRENT_TIMESTEP)
	{
		curts = CURRENT_TIMESTEP;
		tfin = CURRENT_TIME + curts;
	}

	curdelta = tfin - CURRENT_TIME;

	domain = Get_Domain(Oil_Domain);
	gdomain = Get_Domain(Gas_Domain);
	mixdomain = Get_Domain(1);



	thread_loop_c(t, mixdomain)
	{
		begin_c_loop(c, t)
		{
			press_in = C_P(c, t)+101325;
			C_UDMI(c, t, 0) = press_in / 100000;
		}
		end_c_loop(c, t)
	}

	thread_loop_c(t, domain)
	{
		begin_c_loop(c, t)
		{
			C_UDMI(c, t, 3) = C_R(c, t);
			C_UDMI(c, t, 9) = C_MU_L(c, t);
			C_UDMI(c, t, 1) = C_UDSI(c, t, 0);
		}
		end_c_loop(c, t)
	}

	thread_loop_c(gpt, gdomain)
	{
		begin_c_loop(gc, gpt)
		{
			/*C_UDMI(gc, gpt, 1) = C_UDMI(gc, gpt, 2);*/
			press_in = C_UDMI(gc, gpt, 0);
			gvf_in = C_VOF(gc, gpt);
			gio_in = C_UDMI(gc, gpt, 1);
			cv_in = C_VOLUME(gc, gpt);
			odens_in = C_UDMI(gc, gpt, 3);
			ovisc_in = C_UDMI(gc, gpt, 9);
			gdens_in = C_R(gc, gpt);
			gvisc_in = C_MU_L(gc, gpt);
			C_UDMI(gc, gpt, 5) = gdens_in;
			C_UDMI(gc, gpt, 7) = gvf_in;
			datout = changing_pars(press_in, gvf_in, gio_in, cv_in, odens_in, gdens_in, ovisc_in, gvisc_in);
			C_UDMI(gc, gpt, 8) = datout.gvf;
			C_UDMI(gc, gpt, 2) = datout.gio;
			C_UDMI(gc, gpt, 4) = datout.od;
			C_UDMI(gc, gpt, 6) = datout.gd;
			C_UDMI(gc, gpt, 10) = datout.osound;
			C_UDMI(gc, gpt, 11) = datout.gsound;
			C_UDMI(gc, gpt, 12) = datout.ov;
			C_UDMI(gc, gpt, 13) = datout.gv;
		}
		end_c_loop(gc, gpt)
	}

	/*  thread_loop_c(t,domain)
	{
	begin_c_loop (c,t)
	{
	   C_VOF(c, t) = C_UDMI(c, t, 0);
	}
	end_c_loop (c,t)
	}*/

}


DEFINE_ADJUST(adjustment_small, dom)
{	FILE *fil;

	Domain *domain;
	cell_t c;
	Thread *t;
	domain = Get_Domain(Oil_Domain);

	thread_loop_c(t, domain)
	{
		begin_c_loop(c, t)
		{
			C_UDMI(c, t, 1) = C_UDSI(c, t, 0);
		}
		end_c_loop(c, t)
	}

}

DEFINE_LINEARIZED_MASS_TRANSFER(oiltogas_lin, c, t, fpi, fsi, tpi, tsi, fromlt, tolt)
{
	double mtkmax;
	double dt;

	double k1;
	double k2;

	double gd1;
	double gd2;
	double gvf1;
	double gvf2;

	
	if (curdelta<=equil_time)
	{
		dt=equil_time;
	}
	else 
	{
		dt=curdelta;
	}

	gd1 = C_UDMI(c, t, 1);
	gd2 = C_UDMI(c, t, 2);

	gvf1 = C_UDMI(c, t, 7);
	gvf2 = C_UDMI(c, t, 8);

	mtkmax = -C_UDMI(c, t, 17)/dt;
	/*mtkmax = (gd1 - gd2)*gas_surf_density / dt;*/

	if (gd2<gd1)
	{
		k2 = mtkmax / (((1.0 - gvf2)*gvf1) / gvf2 - (1.0 - gvf1));
		k1 = k2*(1.0 - gvf2) / gvf2;
		*fromlt = 0.0;
		*tolt = 0.0;
		/**fromlt = k1 + k2;
		*tolt = -k1 - k2;*/
	}
	else
	{
		*fromlt = 0.0;
		*tolt = 0.0;
		mtkmax = 0.0;
	}
	if (NNULLP(THREAD_STORAGE(t, SV_MT_DS_DP))) C_STORAGE_R(c, t, SV_MT_DS_DP) = 0.0;
	return(mtkmax);
}


DEFINE_LINEARIZED_MASS_TRANSFER(gastooil_lin, c, t, fpi, fsi, tpi, tsi, fromlt, tolt)
{

	double mtkmax;
	double dt;

	double k1;
	double k2;

	double gd1;
	double gd2;
	double gvf1;
	double gvf2;

	if (curdelta<=equil_time)
	{
		dt=equil_time;
	}
	else 
	{
		dt=curdelta;
	}

	gd1=C_UDMI(c,t,1);
	gd2=C_UDMI(c,t,2);

	gvf1=C_UDMI(c,t,7);
	gvf2=C_UDMI(c,t,8);

	mtkmax = C_UDMI(c, t, 17)/dt;
	/*mtkmax = (gd2 - gd1)*gas_surf_density / dt;*/

	if (gd2>gd1)
	{
		k2 = mtkmax / (((1.0 - gvf2)*gvf1) / gvf2 - (1.0 - gvf1));
		k1 = k2*(1.0 - gvf2) / gvf2;
		*fromlt = 0.0;
		*tolt = 0.0;
		/**fromlt = k1 + k2;
		*tolt = -k1 - k2;*/
	}
	else
	{
		*fromlt = 0.0;
		*tolt = 0.0;
        mtkmax = 0.0;
	}
	if(NNULLP(THREAD_STORAGE(t, SV_MT_DS_DP))) C_STORAGE_R(c,t,SV_MT_DS_DP) = 0.0;
	return(mtkmax);
}



DEFINE_MASS_TRANSFER(oiltogas_mtr, c, t, from_phase, from_spec, to_phase, to_spec)
{
	double mtkmax;
	double dt = curdelta;

	double od1;
	double od2;
	double ovf1;
	double ovf2;

	FILE *fil;

	od1=C_UDMI(c,t,3);
	od2=C_UDMI(c,t,4);

	ovf1=1.0-C_UDMI(c,t,7);
	ovf2=1.0-C_UDMI(c,t,8);
	
	if ((curdelta > 0.0)&&(ovf1>ovf2))
	{
		mtkmax = (od1*ovf1 - od2*ovf2) / dt;
	}
	else
	{
		mtkmax = 0.0;
	}

	/*fil = fopen(testouts, "a");
	fprintf(fil, "%7.4f : %7.4f : %7.4f \n", od1, od2, mtkmax);
	fclose(fil);*/

	
    if (mtkmax<0.0) mtkmax=0.0;    
	Message ("value=%f",mtkmax);
	return( mtkmax/50);
}


DEFINE_MASS_TRANSFER(gastooil_mtr, c, t, from_phase, from_spec, to_phase, to_spec)
{
	double mtkmax;
	double dt = curdelta;

	double gd1;
	double gd2;
	double gvf1;
	double gvf2;

	gd1=C_UDMI(c,t,5);
	gd2=C_UDMI(c,t,6);

	gvf1=C_UDMI(c,t,7);
	gvf2=C_UDMI(c,t,8);

	if ((curdelta > 0.0)&&(gvf1>gvf2))
	{
		mtkmax = (gd1*gvf1 - gd2*gvf2) / dt;
	}
	else
	{
		mtkmax = 0.0;
	}

    if (mtkmax<0.0) mtkmax=0.0;

	return(mtkmax);
}


DEFINE_SOURCE(uds_source, c, t, dS, eqn)
{
	real usource;

	real do1;
	real do2;
	real dt;
	real gio1;
	real gio2;

	FILE *fil;

	if (curdelta<=equil_time)
	{
		dt=equil_time;
	}
	else 
	{
		dt=curdelta;
	}

	do1 = C_UDMI(c, t, 3);
	do2 = C_UDMI(c, t, 4);

	gio1 = C_UDMI(c, t, 1);
	gio2 = C_UDMI(c, t, 2);

	if (curdelta > 0.0)
	{
		usource = C_UDMI(c,t,18)/dt;
		/*usource = (gio2 - gio1) / dt;*/
	}
	else
	{
		usource = 0.0;
	}

	/*fil = fopen(testouts, "a");
	fprintf(fil, "%7.4f : %7.4f : %7.4f : %7.4f \n", gio2, gio1, usource, dt);
	fclose(fil);*/

	dS[eqn] = 0.0;
	/*dS[eqn] = (do2 - do1) / dt;*/

	return(usource);
}


DEFINE_EXECUTE_ON_LOADING(loadUDF,libname)
{
    FILE *f;
    char c[2];
	int i = 0;    
	
    /* Reading PVT tables for gas and oil */

    f=fopen(PVT_tables,"r"); 
	fscanf(f,"%le",&oil_surf_density);
    fscanf(f,"%le",&gas_surf_density);
    fscanf(f,"%d",&gas_line_counter);
    for(i=0;i<gas_line_counter;i++) /*gas*/
    {
       fscanf(f,"%le %le %le", &gas_pvt[i][0], &gas_pvt[i][1], &gas_pvt[i][2]);
    }
    fscanf(f,"%s",c);
    while(!feof(f)) /*oil*/
    {
        fscanf(f,"%le %le %le %le", &oil_pvt[oil_line_counter][0], &oil_pvt[oil_line_counter][1], &oil_pvt[oil_line_counter][2], &oil_pvt[oil_line_counter][3]);
        oil_line_counter++;
    }
    fclose(f);
	fillPminmax();
	Message("Loaded...");	
	Message("Loaded...%d", gas_line_counter);
	
}
