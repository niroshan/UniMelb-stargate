
/*
Module Description

This module contains the position control algorithm used for position and 
tracking for Antenna position control systems. This module can be used
as is for integration in a real-time adaptive system or as part of a 
development environment in which the user may change system parameters
interactively and run simulations using a digitised plant model.

Initialisation

To initialise the position controller module the following routines must be
called:

	init_position_controllers ( );

Once initialised, there are other routines that must be called
at certain times to update internal controller variables. These
include:
	reset_pos_parameters ( )
		This resets the PI integrator and loads the measured
		position into controller variables that need to be reset
		to enable smooth changeover.

	reset_pos_filters ( )
		This should be called while slewing to enable smooth changeover
		to closed position loop.

	Remember that driving in rate mode is NOT open-loop but is position
	closed loop and merely tracking a constant ramp trajectory.

Slewing

	When slewing, the user normally applies a constant rate command
	to the drives. This may be done by using the rate_limiter
	routine. By passing in the desired MOTOR speed, not the load speed,
	we can safely accelerate the motor and antenna system. The rate
	limiter is used in BOTH the position loop and slewing loops.
	The rate limiter does BOTH resonance energy filtering and
	rate limiting.

	e.g. SLEW_SPEED = 10500;  / 1750rpm =10500 deg/sec /

			volts = rate_limiter ( SLEW_SPEED );

Positioning

	When in precision position closed loop, one calls the control loop
	routine, passing in the appropriate arguments, and the desired
	MOTOR speed is returned. This is then passed, as above, to the
	rate_limiter.

	e.g.	motor_speed = control_loop ( axis, adapt_switch, cmd, msd );

			volts = rate_limiter ( motor_speed );

	where axis is the axis to control
			adapt_switch turns adaptation on or off
			cmd is the desired commanded position in degrees
			msd is the measured position in degrees

Description of Position Loop

	The position loop is a PI controller with anti-windup running sum
	integrator. There is a second order feed-forward component that
	is used in tracking mode ( when the adaptation switch is on ).

	The first-order derivative feed-forward component is adaptive. It
	is a performance oriented adaptive control whose criterion is to
	reduce position error as well as reduce integrator action in the
	feed-back PI loop. The method used to adapt is the MIT-Rule
	which is merely a gradient type estimator. There are a number of
	parameters that must be setup for this and for the PI loop.

	All of the parameters that need to be setup may be found in
	the structure par[axis].

Structures

	There are 4 main structures that allow one to look into the controller
	operation. This is particularly handy when debugging as all variables
        may be inspected in structure con [axis].

	controller_parameter_struct par [ 2 ];
		Contains controller parameters. e.g. proportional and integral
		gains, adaptive controller constants. All are constant
		until changed in NVRAM. Note

	controller_struct con [ 2 ];
		Contains all temporary variables that belong to each axis.

	res_filter_struct res_filter [ 2 ];
		Contains regressor and coefficient values for the anti-resonance
		filters.

	diff_filter_struct diff_filter [ 2 ];
		Contains the regressor and coefficient values for the
		differential filters used in the feed-forward.

	mit_adaptation_struct mit [ 2 ];
		Contains all temporary variables used in the adaptation routine.

*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "poscon.h"

#include "seq.h"
#include "utilities.h"
#include "general.h"

#define DAMP 0.6
#define pi 3.141592654


/*---- Azimuth Prototypes -----*/
static void init_position_controller (unsigned int axis);
static void ff_adaptation (unsigned int axis, double cmd);
static void pid_controller (unsigned int axis,
                            unsigned int adaptation_switch,
                            double cmd, double msd );
static double rate_to_volts (unsigned int axis, double az_control);
static int get_controller_parameters (FILE *fc, unsigned int axis);

/*
static double resonance_filter ( unsigned int axis, double filter_in );
*/


/*-------------------- Local Variables for Module ------------------------*/

/*---- General Variables ----*/

		/*---- Structures (typedefed in posext.h ) ----*/

/* Control algorithm constants read in from NVRAM */
static controller_parameter_struct par [2];

/* Structures that allow generic controller routines */

static controller_struct con [2];
static res_filter_struct res_filter [2];
static diff_filter_struct diff_filter [2];

static mit_adaptation_struct mit [2];


/*--------------------------------- Module Start -------------------------*/

/***************************************************************************
**
** init_position_controllers ().
**
**	Routine to initialise both axes position controllers.
**
** Reads the values from the file "controller.data".
**
** Parameters:
**	None
**
***************************************************************************/

int init_position_controllers () 
{
   int status1, status2;

   FILE *fc;

   #ifdef DEBUG
   debug_going_down ();
   print_spaces ();
   puts ("In init_position_controllers()");
   #endif

   if ((fc = fopen ("controller.data", "r")) == NULL) {
      printf ("\nFatal Error: cannot access the file containing "
            "initial controller values.\n\n");
      exit (EXIT_FAILURE);
   }


   /* Get controller Parameters from the file "controller.data" */
   status1 = get_controller_parameters (fc, AZ_AXIS);
   status2 = get_controller_parameters (fc, EL_AXIS);

   fclose (fc);

   init_position_controller (AZ_AXIS);
   init_position_controller (EL_AXIS);

   #ifdef DEBUG
   debug_going_up ();
   #endif

   return (status1 | status2);

} /* end of init_position_controllers */


/***************************************************************************
**
** reset_pos_parameters ().
** 
**      Routine to reset the position controller parameters for azimuth 
** drives. This is necessary to ensure smooth changeover to precision 
** positioning. The regressor vectors ( previous values ) in the controller
** filters have to reset to avoid transients.
**
** Parameters:
**	el - Pointer to axis structure for elevation
**
************************************************************************* */

void reset_pos_parameters (unsigned int axis, double msd)
{
   /*---- Reset some controller Parameters */
   con [axis].cmd_prev = msd;
   con [axis].cmd_mit_prev = msd;

   con [axis].rate_volts = 0;
   con [axis].cmd_diff_1 = 0;
   con [axis].cmd_dotdot_1 = 0;
   con [axis].cmd_dot = 0;

   /* pid_controller */
   con [axis].int_sum = 0;
   con [axis].control_limited = 0;

   /* Slew rate limiter */
   con [axis].dem_rate_1 = 0;

   /* Initialise the ff gain */
   mit [axis].ff_gain = par [axis].initial_mit_ff_gain  ;

   /* Set Resonance filter regressors */
   res_filter [axis].ym_2 = 0;
   res_filter [axis].ym_1 = 0;
   res_filter [axis].um_1 = 0;

} /* end of reset_pos_parameters */


/**************************************************************************
**
** reset_pos_filters ().
**
**	Routine to reset the position controller parameters for azimuth
** drives. This is necessary to ensure smooth changeover to precision
** positioning. The regressor vectors (previous values) in the controller
** filters have to reset to avoid transients.
**
**      This routine should be called when slewing to ensure that when we 
** change over to position loop mode we have smooth changeover.
**
** Parameters:
** 	axis
**	el - msd position
**
** Returns:
**	None
**
**************************************************************************/

void reset_pos_filters (unsigned int axis, double msd)
{
   /*---- Reset some controller Parameters */
   con [axis].cmd_prev = msd;
   con [axis].cmd_mit_prev = msd;

   con [axis].cmd_diff_1 = 0;
   con [axis].cmd_dotdot_1 = 0;
   con [axis].cmd_dot = 0;

   /* pid_controller */
   con [axis].int_sum = 0;

   /* Initialise the ff gain */
   mit [axis].ff_gain = par [axis].initial_mit_ff_gain  ;

} /* end of reset_pos_filters */


/*************************************************************************
** init_position_controller ().
** 
**      Routine to initialise local variables in this module. This routine
** should be called prior to carrying out real-time control and ideally 
** before applying power to the drive system. This may be an important 
** issue when attempting to avoid the situation where we may have initialised
** the last position value, gone into slew mode and then gone into
** position mode.
**
**      Because this routine modifies controller filter regressors, it is
** essential to call this routine after a demanded position has been 
** determined.
**
** Parameters:
**	axis - axis number AZ_AXIS or EL_AXIS
**
** Returns:
**	None
**
*************************************************************************/

static void init_position_controller (unsigned int axis) 
{
   double a,w,b,g,wr;

   /* for resonance filter */
   con [axis].control_1 = 0;

   /* pid_controller */
   con [axis].int_sum = 0;
   con [axis].cmd_prev = 0;
   con [axis].cmd_dot = 0;
   con [axis].control_limited = 0;

   /* Slew rate limiter */
   con [axis].control_rate_1 = 0;
   con [axis].dem_rate_1 = 0;

   /* Get Resonance filter coefficients */

   /*
   ** This WAS the first order approximation. I left it here in case
   ** we had trouble in the field.
   */
   /*
   res_filter [axis].alpha 
             = exp (-par [axis].res_bw * CONTROL_PERIOD * 2 * pi);
   res_filter [axis].beta = (1 - res_filter [axis].alpha);
   */

   /* Resonance ZOH Equivalent */
   /* Second order filter */
   wr = par [axis].res_bw * 2.0 * pi;
   a = exp (-DAMP * wr * CONTROL_PERIOD);
   w = wr * sqrt (1 - DAMP * DAMP);
   b = cos (w * CONTROL_PERIOD);
   g = sin (w * CONTROL_PERIOD);

   res_filter [axis].b1 = -(1 - a * (b + (DAMP * wr * g / w)));
   res_filter [axis].b2 = -(a * a + a * ((DAMP * wr * g / w) - b));
   res_filter [axis].a1 = -2 * a * b;
   res_filter [axis].a2 = a * a;

   /* Get differentiator Filter coefficients */
   diff_filter [axis].alpha 
                = exp (-par [axis].diff_bw * CONTROL_PERIOD * 2 * pi );
   diff_filter [axis].beta = (1 - diff_filter [axis].alpha);

} /* end of init_position_controller */


/**************************************************************************
** 
**  pid_control_loop ().
**  
**	Routine that does does the control loop. It calculates the measured
** error decides on whether to do precision control or rate control. If it
** does precision control then it calls the PID and adaptive algorithms.
**
**	This routine should be called by the external module that wishes to 
** do closed-loop control.
**
**	It is possible to call the precision position loop with adaptation
** turned off. This may be required in certain modes that want to do small 
** step perturbations.
**
** Parameters:
**	(unsigned int)axis = axis to be controlled
**	(unsigned int)adaptation_switch : = 0 turn adaptation off and
**	        do not use feed-forward.
**	(double)cmd = commanded az position in degrees
**	(double)msd = measured az position in degrees
**
** Returns:
**	(double) control deg/sec
**
**************************************************************************/

double pid_control_loop (unsigned int axis,
			 unsigned int adaptation_switch,
			 double cmd, double msd)
{
   /*---- Measured Position Error */
   con [axis].err = cmd - msd;

   /*---- Do precision loop ----*/
   pid_controller (axis, adaptation_switch, cmd, msd);

   if (adaptation_switch)
      ff_adaptation (axis, cmd);

   return (con [axis].control_unfilt);

} /* end of pid_control_loop */


/**************************************************************************
**
** ff_adaptation ().
**
** 	This routine carries out the MIT-Rule adaptation as described 
** in the references. Briefly, the aim of the algorithm is to dynamically
** adjust the position feed-forward parameter to reduce measured Antenna
** error with respect to the demanded position. It does this by integrating
** the sum of the position error and the integral error. By using the 
** integral error signal derived from the PID controller, it is possible
** to reduce the energy stored in the integral action.
**
**	Tuning the adaptive algorithm is particularly important and can have
** a marked affect on the performance of the Antenna. Careless choice of the
** 'alpha' and 'beta' coefficients may in fact result in instability or
** oscillations. Some guidelines will be presented in late 1990 as to the
** tuning procedure.
**
** Parameter:
**	axis = axis to control
**	cmd = commanded position
**
**************************************************************************/

static void ff_adaptation (unsigned int axis, double az_cmd)
{
   /* Do differentiator */
   con [axis].cmd_dot = (az_cmd - mit [axis].cmd_prev) / CONTROL_PERIOD;
   mit [axis].cmd_prev = az_cmd;

   /* Normalise for one rate */
   mit [axis].norm = (par [axis].mit_maxref) 
                / (fabs (con [axis].cmd_dot) + 1e-10);

   mit [axis].sign_cmd_dot = (con [axis].cmd_dot 
                / (fabs (con [axis].cmd_dot) + 1e-10));

   mit [axis].alpha_err = par [axis].mit_alpha 
                * mit [axis].norm * con [axis].err;

   mit [axis].beta_err = par [axis].mit_beta * con [axis].int_sum;

   mit [axis].incr = -(mit [axis].alpha_err + mit [axis].beta_err)
                * mit [ axis ].sign_cmd_dot;

   /* Apply some adaptation limiting for robustness */
   if (mit [axis].incr > par [axis].mit_maxincr)
      mit [axis].incr = par [axis].mit_maxincr;

   if (mit [axis].incr < -par [axis].mit_maxincr)
      mit [axis].incr	= -par [axis].mit_maxincr;

   mit [axis].ff_gain = mit [axis].ff_gain - mit [axis].incr;

   /* Apply limits to the adapted feed-forward parameter */
   if (mit [axis].ff_gain > par [axis].mit_ff_high)
      mit [axis].ff_gain = par [axis].mit_ff_high;
   if ( mit [axis].ff_gain < par [axis].mit_ff_low)
      mit [axis].ff_gain = par [axis].mit_ff_low;

} /* end of mit_rule_adaptation */



/**************************************************************************
**
** pid_controller ().
** 
**      Routine to carry out PID feed-back and Feed-Forward Control. 
** Although the original version did not incorporate derivative action, 
** it has been found that some D action may help to reduce the effects of
** stiction.
**
** Calls:
**	resonance_filter
**
** Arguments:
**	axis
**	adaptation switch
**	cmd
**	msd
**
*************************************************************************/

static void pid_controller (unsigned int axis,
			    unsigned int adaptation_switch,
			    double cmd,	double msd) 
{
   /*---- Differentiate the commanded position */
   con [axis].cmd_diff = (cmd - con [axis].cmd_prev) / CONTROL_PERIOD;
   con [axis].cmd_prev = cmd;

   if ( con [axis].cmd_diff > 0.3 )
      con [axis].cmd_diff = 0.3;
   if ( con [axis].cmd_diff < -0.3 )
      con [axis].cmd_diff = -0.3;

   con [axis].cmd_dotdot 
      = (con [axis].cmd_diff - con [axis].cmd_dotdot_1) / CONTROL_PERIOD;

   con [axis].cmd_dotdot_1 = con [axis].cmd_diff;

   if (con [axis].cmd_dotdot > 0.2)
      con [axis].cmd_dotdot = 0.2;
   if (con [axis].cmd_dotdot < -0.2)
      con [axis].cmd_dotdot = -0.2;

   con [axis].cmd_diff = diff_filter [axis].beta * con [axis].cmd_diff 
      + con [axis].cmd_diff_1 * diff_filter [axis].alpha;
   con [axis].cmd_diff_1 = con [axis].cmd_diff;

   if (fabs (con [axis].err) < 0.2)
      con [axis].err = 0;

   /*---- Do integrator */
   con [axis].int_sum = con [axis].int_sum + con [axis].err;

   if (con [axis].int_sum > par [axis].int_positive_limit)
      con [axis].int_sum = par [axis].int_positive_limit;
   if (con [axis].int_sum < par [axis].int_negative_limit)
      con [axis].int_sum = par [axis].int_negative_limit;

   con [axis].int_action = par [axis].integral_gain * con [axis].int_sum;

   /*---- Do feed-forward on reference and find new error */
   if (adaptation_switch == 0)
      con [axis].pred = cmd;
   else
      con [axis].pred = cmd + mit [axis].ff_gain * con [axis].cmd_diff 
        + par [axis].parabolic_gain * con [axis].cmd_dotdot;

   con [axis].control_unfilt = par [axis].prop_gain * (con [axis].pred - msd) 
        + con [axis].int_action;

} /* end of pid_controller */


/**************************************************************************
**
** resonance_filter ().
**
**      Routine to filter control signal to avoid exciting resonance in the
** Antenna structure.
**
** 1st or 2nd order.
**
**************************************************************************/

/*
double resonance_filter (unsigned int axis, double filter_in) 
{
   double control_tmp;

   control_tmp = (-res_filter [axis].a2 * res_filter [axis].ym_2 
        - res_filter [axis].a1 * res_filter [axis].ym_1 
        - res_filter [axis].b1 * filter_in 
        + res_filter [axis].b2 * res_filter [axis].um_1);

   res_filter [axis].ym_2 = res_filter [axis].ym_1;
   res_filter [axis].ym_1 = control_tmp;
   res_filter [axis].um_1 = -filter_in;

   return (control_tmp);
} 
*/


/***************************************************************************
**
** rate_limiter ().
**
**      Routine to limit rate action applied to the motor drives. This is
** primarily designed to avoid damaging the motor/gear system. Also scales
** and outputs the control to the DACs.
**
**      This routine also calls the resonance filter before applying the
** rate limits. Note that the limiting was done after the filter in case
** the filter caused large spikes to enter the system.
**
**      The maximum acceleration is specified in terms of load speed (this
** could be made motor side if desired), and if we exceed this limit over
** the control period, we increment theprevious control value (.dem_rate_1)
** by the maximum allowable rate in the control period (that is whay we
** multiply by the CONTROL_PERIOD).
**
**      Also we reset the integrator when we limit.
**      
** Parameters:
**	double dem_rate - demanded rate in deg/sec
**
** Returns:
**	(double)voltage to be applied to drives
**
***************************************************************************/

double rate_limiter (unsigned int axis, double dem_rate)
{
   double control, volts; /* max_acc; */

   /*max_acc = par [axis].max_acceleration * par [axis].gear_ratio; */

   /*---- Smooth demanded rate using a second order filter to ensure
     we do not excite resonances. Removed the call to this
     routine on 4/4/95 (ST and SRW)
   */

   /* control = resonance_filter (axis, dem_rate); */
   control = dem_rate;

   /* Only apply rate limit outside a small rate region */
   /*
   if (fabs (control) >= par [axis].rate_limit_region) {
      if ((control - con [axis].dem_rate_1) / CONTROL_PERIOD > max_acc) {
         control = con [axis].dem_rate_1 + max_acc * CONTROL_PERIOD;
         con [axis].int_sum = 0;
      }
      else
         if ((control - con [axis].dem_rate_1) / CONTROL_PERIOD < -max_acc) {
            control = con [axis].dem_rate_1 - max_acc * CONTROL_PERIOD;
            con [axis].int_sum = 0;
         }
   }
   */

   /*---- This value is used by the control output routines */
   con [axis].control_limited = control;

   /*---- Form Voltage and clip if necessary. N.B. This routine
     does NOT do I/O. It merely determines voltage. */
   volts = rate_to_volts (axis, control);
   con [axis].dem_rate_1 = volts * par [axis].rate_gain * par [axis].sign;

   return (volts);

} /* end of rate_limiter */


/***************************************************************************
**
** rate_to_volts ().
**
**      Routine to convert rate to volts.
**      
**      This routine merely forms the DAC word to be output at next control
** period.
**
** parameters:
**	axis    - az or el
**	control - deg/sec
**
** Returns:
**	(double)volts to be applied to drives
**
***************************************************************************/

static double rate_to_volts ( unsigned int axis, double control_out )
{

   /* we should output the last value to the Azimuth DAC */
   con [axis].rate_volts = (control_out / par [axis].rate_gain);

   if (con [axis].rate_volts > par [axis].max_volts)
      con [axis].rate_volts = par [axis].max_volts;
   if (con [axis].rate_volts < par [axis].min_volts)
      con [axis].rate_volts = par [axis].min_volts;

   /* Apply drive direction sign */
   con [axis].rate_volts = par [axis].sign * con [axis].rate_volts;

   return (con [axis].rate_volts);

} /* end of rate_to_volts */


/*************************************************************************
** get_controller_parameters ().
**
** Sub-function of init_position_controller ().
**
*************************************************************************/

int get_controller_parameters (FILE *fc, unsigned int axis)
{
   char line[LINE];

   fgets (line,LINE, fc);
   sscanf (line, "%lf", &par[axis].prop_gain);

   fgets (line,LINE, fc);
   sscanf (line, "%lf", &par[axis].integral_gain);

   fgets (line,LINE, fc);
   sscanf (line, "%lf", &par[axis].res_bw);

   fgets (line,LINE, fc);
   sscanf (line, "%lf", &par[axis].diff_bw);

   fgets (line,LINE, fc); /* line 5 */
   sscanf (line, "%lf", &par[axis].int_positive_limit);

   fgets (line,LINE, fc);
   sscanf (line, "%lf", &par[axis].int_negative_limit);

   fgets (line,LINE, fc);
   sscanf (line, "%d", &par[axis].sign);

   fgets (line,LINE, fc);
   sscanf (line, "%lf", &par[axis].mit_alpha);

   fgets (line,LINE, fc);
   sscanf (line, "%lf", &par[axis].mit_beta);

   fgets (line,LINE, fc); /* line 10 */
   sscanf (line, "%lf", &par[axis].mit_maxref);

   fgets (line,LINE, fc);
   sscanf (line, "%lf", &par[axis].mit_maxincr);

   fgets (line,LINE, fc);
   sscanf (line, "%lf", &par[axis].mit_ff_high);

   fgets (line,LINE, fc);
   sscanf (line, "%lf", &par[axis].mit_ff_low);

   fgets (line,LINE, fc);
   sscanf (line, "%lf", &par[axis].initial_mit_ff_gain);

   fgets (line,LINE, fc); /* line 15 */
   sscanf (line, "%lf", &par[axis].parabolic_gain);

   fgets (line,LINE, fc);
   sscanf (line, "%lf", &par[axis].max_volts);

   fgets (line,LINE, fc);
   sscanf (line, "%lf", &par[axis].min_volts);

   fgets (line,LINE, fc);
   sscanf (line, "%lf", &par[axis].max_acceleration);

   fgets (line,LINE, fc);
   sscanf (line, "%lf", &par[axis].gear_ratio);

   fgets (line,LINE, fc); /* line 20 */
   sscanf (line, "%lf", &par[axis].rate_gain);

   fgets (line,LINE, fc);
   sscanf (line, "%lf", &par[axis].rate_limit_region);

   fgets (line,LINE, fc);
   sscanf (line, "%lf", &par[axis].slew_region);

   return 1;

} /* end of get_controller_parameters */


double get_slew_region ()
{
   return par [0].slew_region;
}
