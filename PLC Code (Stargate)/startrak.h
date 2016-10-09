extern void star_track (double star_ra, double star_dec, int epoch,
			int day, int month, int year, int hour, int minute,
			double second,
			double station_long, double station_lat,
			double *azcom, double *elcom);

extern double juldat (int day, int month, int year,
			  int hour, int minute, double second);
