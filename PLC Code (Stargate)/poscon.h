/***************************************************************************

 POSCON.H

***************************************************************************/

#ifndef __POSCON_H
#define __POSCON_H


#define AZ_AXIS 0
#define EL_AXIS 1


/* Controller parameters structure. Holds the controller constants.
   and is usually read from from a file. For both axes.
   approximately 22 values per axis.
*/

typedef struct {
   double prop_gain;		/* Controller Proportional gain (Motor deg)*/
   double integral_gain;	/* Integral Gain */
   double res_bw;		/* Resonance Filter cut-off frequency Hz */
   double diff_bw;		/* Feed-forward differentiator filter cut-off */
   double int_positive_limit;	/* Integrator positive limit */
   double int_negative_limit;	/* Integrator negative limit */
   int sign;			/* Motor drive sign. +V*sign -> rate */
   double mit_alpha;		/* MIT-Rule alpha gain */
   double mit_beta;		/* MIT-Rule beta gain */
   double mit_maxref;		/* MIT-Rule Maximum tracking reference */
   double mit_maxincr;		/* Maximum allowed adaptation increment */
   double mit_ff_high;		/* Feed-forward adaptation limit High */
   double mit_ff_low;		/* Feed-forward adaptation limit Low */
   double initial_mit_ff_gain;	/* Initial Feed-forward gain. Reset value */
   double parabolic_gain;	/* Second order feed-forward coefficient */
   double max_volts;		/* Full-speed volts on drive */
   double min_volts;		/* -Full speed volts (probably could lump 
                                   with max_volts */
   double max_acceleration;	/* maximum acceleration of drive in deg/s/s */
   double gear_ratio;		/* Gear ratio N */
   double rate_gain;		/* Motor drive rate gain. deg/sec/V */
   double rate_limit_region;    /* region where rate limiting begins */
   double slew_region;		/* used by sequencer to determine slew region */

} controller_parameter_struct;


/* Controller variables structure */
typedef struct {
   double err;			/* measured - commanded */
   double cmd_dot;
   double pred;
   double control_unfilt;
   double dem_rate;
   double dem_rate_1;
   double cmd_prev;	        /* previous reference command */
   double cmd_mit_prev;
   double rate_volts;
   double cmd_diff_1;
   double cmd_dotdot_1;
   double int_sum;		/* integrator summation */
   double control_limited;
   double control_rate_1;
   double control_1;
   double cmd_diff;
   double cmd_dotdot;
   double int_action;		/* integrator action */

   int rate_dac;		/* DAC word */
} controller_struct;


/* Anti-resonance filter structure. Holds regressor and coefficient values. */
typedef struct {
   double ym_2;
   double ym_1;
   double um_1;
   double a1;
   double a2;
   double b1;
   double b2;
} res_filter_struct;


/* Command differentiator filter structure */
typedef struct {
   double alpha;
   double beta;
} diff_filter_struct;


/* MIT-Rule structure */
typedef struct {
   double norm;
   double alpha_err;
   double beta_err;
   double incr;
   double ff_gain;
   double cmd_prev;
   double sign_cmd_dot;
} mit_adaptation_struct;


/*
  This routine should be called by the sequencer
  when sitting in standby to allow smooth changeover.
*/
extern void reset_pos_parameters (unsigned int axis, double msd);

/*
  Reset filter coefficients only in slewing mode to allow
  smooth changeover to position loop.
*/
extern void reset_pos_filters (unsigned int axis, double msd);

/*
  Routine to initialise position controllers. This should be done
  at system startup as well as after changing controller
  parameters in the file.
*/
extern int init_position_controllers (void);

/* Control Position loop. Returns the MOTOR demanded speed */
extern double pid_control_loop (unsigned int axis,
      unsigned int adaptation_switch,
      double cmd, double msd);

/*
  rate and resonance filter. Should be called in BOTH
  position and slewing modes. dem_rate is MOTOR speed
  in deg/sec. e.g. 1750rpm = 10500 deg/sec.
  Returns the voltage required to be placed on the drives.
*/
extern double rate_limiter (unsigned int axis, double dem_rate);

double get_slew_region (void);

#endif

/*---------------------------- End of poscon.h ---------------------------*/
