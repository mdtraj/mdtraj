#ifndef _trr_header_h_
#define _trr_header_h_

typedef struct		/* This struct describes the order and the	*/
/* sizes of the structs in a trjfile, sizes are given in bytes.	*/
{
    int  bDouble;        /* Double precision?                            */
    int	ir_size;	/* Backward compatibility		        */
    int	e_size;		/* Backward compatibility		        */
    int	box_size;	/* Non zero if a box is present			*/
    int   vir_size;       /* Backward compatibility		        */
    int   pres_size;      /* Backward compatibility		        */
    int	top_size;	/* Backward compatibility		        */
    int	sym_size;	/* Backward compatibility		        */
    int	x_size;		/* Non zero if coordinates are present		*/
    int	v_size;		/* Non zero if velocities are present		*/
    int	f_size;		/* Non zero if forces are present		*/

    int	natoms;		/* The total number of atoms			*/
    int	step;		/* Current step number				*/
    int	nre;		/* Backward compatibility		        */
    float	tf;		/* Current time					*/
    float	lambdaf;		/* Current value of lambda			*/
    double	td;		/* Current time					*/
    double	lambdad;		/* Current value of lambda			*/
} t_trnheader;

int do_trnheader(XDRFILE *xd, char bRead, t_trnheader *sh);

#endif
