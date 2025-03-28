#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_math.h>

/* Enums will be used for indexing purposes.   */
enum stateVariable { _RiiP, _RiiP_cAMP, _RiiP_C, _RiiP_C_cAMP, _C, _Rii_cAMP, _Rii_C_cAMP, _RiiP_CaN, _RiiP_cAMP_CaN, _AKAR4_C, _AKAR4p, numStateVar };
enum param { _kf_Rii_C__RiiP_C, _kf_RiiP_CxcAMP__RiiP_C_cAMP, _kf_RiiP_cAMPxC__RiiP_C_cAMP, _kb_RiiP_cAMPxC__RiiP_C_cAMP, _kb_RiiPXcAMP__RiiP_cAMP, _kf_RiiPXcAMP__RiiP_cAMP, _kf_RiiPxC__RiiP_C, _kb_RiiPxC__RiiP_C, _kf_cAMPxRii__Rii_cAMP, _kb_cAMPxRii__Rii_cAMP, _kf_Rii_CxcAMP__Rii_C_cAMP, _kb_Rii_CxcAMP__Rii_C_cAMP, _kf_RiixC__Rii_C, _kf_Rii_cAMPxC__Rii_C_cAMP, _kb_Rii_cAMPxC__Rii_C_cAMP, _kf_Rii_C_cAMP__RiiP_C_cAMP, _kb_RiixC__Rii_C, _AKAPoff_1, _AKAPoff_3, _AKAPon_1, _AKAPon_3, _kf_C_AKAR4, _kb_C_AKAR4, _kcat_AKARp, _kmOFF, _kmON, _KD_T, _b_AKAP, _AKAR4_ConservedConst, _CaN_ConservedConst, _Rii_C_ConservedConst, _cAMP_ConservedConst, _Rii_ConservedConst, numParam };
enum func { _AKAR4pOUT, numFunc };

/* The error codes indicate how many values a function returns.                             */
/* Each function expects the output buffer to be allocated with at least that many values   */

/* ODE vector field: y' = f(t,y;p)   */
int AKAP79_vf(double t, const double y_[], double *f_, void *par){
	double *p_=par;
	if (!y_ || !f_) return 11;
/* 	constants   */
/* 	parameter values   */
	double kf_Rii_C__RiiP_C = p_[_kf_Rii_C__RiiP_C];                    /* [  0] */
	double kf_RiiP_CxcAMP__RiiP_C_cAMP = p_[_kf_RiiP_CxcAMP__RiiP_C_cAMP];       /* [  1] */
	double kf_RiiP_cAMPxC__RiiP_C_cAMP = p_[_kf_RiiP_cAMPxC__RiiP_C_cAMP];       /* [  2] */
	double kb_RiiP_cAMPxC__RiiP_C_cAMP = p_[_kb_RiiP_cAMPxC__RiiP_C_cAMP];       /* [  3] */
	double kb_RiiPXcAMP__RiiP_cAMP = p_[_kb_RiiPXcAMP__RiiP_cAMP];       /* [  4] */
	double kf_RiiPXcAMP__RiiP_cAMP = p_[_kf_RiiPXcAMP__RiiP_cAMP];       /* [  5] */
	double kf_RiiPxC__RiiP_C = p_[_kf_RiiPxC__RiiP_C];                  /* [  6] */
	double kb_RiiPxC__RiiP_C = p_[_kb_RiiPxC__RiiP_C];                  /* [  7] */
	double kf_cAMPxRii__Rii_cAMP = p_[_kf_cAMPxRii__Rii_cAMP];          /* [  8] */
	double kb_cAMPxRii__Rii_cAMP = p_[_kb_cAMPxRii__Rii_cAMP];          /* [  9] */
	double kf_Rii_CxcAMP__Rii_C_cAMP = p_[_kf_Rii_CxcAMP__Rii_C_cAMP];       /* [ 10] */
	double kb_Rii_CxcAMP__Rii_C_cAMP = p_[_kb_Rii_CxcAMP__Rii_C_cAMP];       /* [ 11] */
	double kf_RiixC__Rii_C = p_[_kf_RiixC__Rii_C];                      /* [ 12] */
	double kf_Rii_cAMPxC__Rii_C_cAMP = p_[_kf_Rii_cAMPxC__Rii_C_cAMP];       /* [ 13] */
	double kb_Rii_cAMPxC__Rii_C_cAMP = p_[_kb_Rii_cAMPxC__Rii_C_cAMP];       /* [ 14] */
	double kf_Rii_C_cAMP__RiiP_C_cAMP = p_[_kf_Rii_C_cAMP__RiiP_C_cAMP];       /* [ 15] */
	double kb_RiixC__Rii_C = p_[_kb_RiixC__Rii_C];                      /* [ 16] */
	double AKAPoff_1 = p_[_AKAPoff_1];                                  /* [ 17] */
	double AKAPoff_3 = p_[_AKAPoff_3];                                  /* [ 18] */
	double AKAPon_1 = p_[_AKAPon_1];                                    /* [ 19] */
	double AKAPon_3 = p_[_AKAPon_3];                                    /* [ 20] */
	double kf_C_AKAR4 = p_[_kf_C_AKAR4];                                /* [ 21] */
	double kb_C_AKAR4 = p_[_kb_C_AKAR4];                                /* [ 22] */
	double kcat_AKARp = p_[_kcat_AKARp];                                /* [ 23] */
	double kmOFF = p_[_kmOFF];                                          /* [ 24] */
	double kmON = p_[_kmON];                                            /* [ 25] */
	double KD_T = p_[_KD_T];                                            /* [ 26] */
	double b_AKAP = p_[_b_AKAP];                                        /* [ 27] */
	double AKAR4_ConservedConst = p_[_AKAR4_ConservedConst];            /* [ 28] */
	double CaN_ConservedConst = p_[_CaN_ConservedConst];                /* [ 29] */
	double Rii_C_ConservedConst = p_[_Rii_C_ConservedConst];            /* [ 30] */
	double cAMP_ConservedConst = p_[_cAMP_ConservedConst];              /* [ 31] */
	double Rii_ConservedConst = p_[_Rii_ConservedConst];                /* [ 32] */
/* 	state variables   */
	double RiiP = y_[_RiiP];                                            /* [  0] */
	double RiiP_cAMP = y_[_RiiP_cAMP];                                  /* [  1] */
	double RiiP_C = y_[_RiiP_C];                                        /* [  2] */
	double RiiP_C_cAMP = y_[_RiiP_C_cAMP];                              /* [  3] */
	double C = y_[_C];                                                  /* [  4] */
	double Rii_cAMP = y_[_Rii_cAMP];                                    /* [  5] */
	double Rii_C_cAMP = y_[_Rii_C_cAMP];                                /* [  6] */
	double RiiP_CaN = y_[_RiiP_CaN];                                    /* [  7] */
	double RiiP_cAMP_CaN = y_[_RiiP_cAMP_CaN];                          /* [  8] */
	double AKAR4_C = y_[_AKAR4_C];                                      /* [  9] */
	double AKAR4p = y_[_AKAR4p];                                        /* [ 10] */
/* 	expressions   */
	double kf_RiiP_cAMP_CaN__CaNXRii_cAMP = b_AKAP * AKAPon_1 + (1 - b_AKAP) * AKAPoff_1;
	double kb_RiiPxCaN__RiiP_CaN = b_AKAP*AKAPon_3  +  (1 - b_AKAP)* AKAPoff_3;
	double kf_RiiP_CaN__RiixCaN = b_AKAP * AKAPon_1 + (1 - b_AKAP) * AKAPoff_1;
	double kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN = b_AKAP*AKAPon_3  +  (1 - b_AKAP)* AKAPoff_3;
	double kf_RiiPxCaN__RiiP_CaN = b_AKAP * ((kb_RiiPxCaN__RiiP_CaN + kf_RiiP_cAMP_CaN__CaNXRii_cAMP)/kmON ) + (1 - b_AKAP) * (kb_RiiPxCaN__RiiP_CaN + kf_RiiP_cAMP_CaN__CaNXRii_cAMP) / kmOFF;
	double kf_CaNxRiiP_cAMP__RiiP_cAMP_CaN = b_AKAP * ((kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN + kf_RiiP_CaN__RiixCaN)/kmON) + (1 - b_AKAP) * (kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN + kf_RiiP_CaN__RiixCaN)/kmOFF;
	double kb_RiiP_CxcAMP__RiiP_C_cAMP = kf_RiiP_CxcAMP__RiiP_C_cAMP * KD_T;
	double AKAR4 = (AKAR4_ConservedConst - (AKAR4_C+AKAR4p));
	double CaN = (CaN_ConservedConst - (RiiP_CaN+RiiP_cAMP_CaN));
	double Rii_C = (Rii_C_ConservedConst - (RiiP_C+RiiP_C_cAMP+C+Rii_C_cAMP+AKAR4_C));
	double cAMP = (cAMP_ConservedConst - (RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN));
	double Rii = (Rii_ConservedConst - (RiiP+RiiP_cAMP+Rii_cAMP+RiiP_CaN+RiiP_cAMP_CaN-C-AKAR4_C));
	double reaction_51 = kf_Rii_C__RiiP_C*Rii_C;
	double reaction_14 = kf_RiiPxC__RiiP_C*RiiP*C - kb_RiiPxC__RiiP_C*RiiP_C;
	double reaction_12 = kf_RiiP_CxcAMP__RiiP_C_cAMP*RiiP_C*cAMP - kb_RiiP_CxcAMP__RiiP_C_cAMP*RiiP_C_cAMP;
	double reaction_43 = kf_RiiPXcAMP__RiiP_cAMP*cAMP*RiiP - kb_RiiPXcAMP__RiiP_cAMP*RiiP_cAMP;
	double reaction_23 = kf_RiiP_cAMPxC__RiiP_C_cAMP*RiiP_cAMP*C - kb_RiiP_cAMPxC__RiiP_C_cAMP*RiiP_C_cAMP;
	double reaction_78 = kf_cAMPxRii__Rii_cAMP*cAMP*Rii - kb_cAMPxRii__Rii_cAMP*Rii_cAMP;
	double reaction_56 = kf_Rii_CxcAMP__Rii_C_cAMP*Rii_C*cAMP - kb_Rii_CxcAMP__Rii_C_cAMP*Rii_C_cAMP;
	double reaction_76 = kf_Rii_cAMPxC__Rii_C_cAMP*Rii_cAMP*C - kb_Rii_cAMPxC__Rii_C_cAMP*Rii_C_cAMP;
	double reaction_62 = kf_Rii_C_cAMP__RiiP_C_cAMP*Rii_C_cAMP;
	double reaction_58 = kf_RiixC__Rii_C*Rii*C - kb_RiixC__Rii_C*Rii_C;
	double reaction_44 = kf_RiiPxCaN__RiiP_CaN*RiiP*CaN - kb_RiiPxCaN__RiiP_CaN*RiiP_CaN;
	double reaction_33 = kf_CaNxRiiP_cAMP__RiiP_cAMP_CaN*CaN*RiiP_cAMP - kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN*RiiP_cAMP_CaN;
	double reaction_48 = kf_RiiP_CaN__RiixCaN*RiiP_CaN;
	double reaction_37 = kf_RiiP_cAMP_CaN__CaNXRii_cAMP*RiiP_cAMP_CaN;
	double reaction_1 = kf_C_AKAR4*C*AKAR4 - kb_C_AKAR4*AKAR4_C;
	double reaction_2 = kcat_AKARp*AKAR4_C;
	memset(f_,0,sizeof(double)*11); /* initialize with 0.0 */
	f_[0] = -reaction_14-reaction_43-reaction_44;
	f_[1] = +reaction_43-reaction_23-reaction_33;
	f_[2] = +reaction_51+reaction_14-reaction_12;
	f_[3] = +reaction_12+reaction_23+reaction_62;
	f_[4] = -reaction_14-reaction_23-reaction_76-reaction_58-reaction_1+reaction_2;
	f_[5] = +reaction_78-reaction_76+reaction_37;
	f_[6] = +reaction_56+reaction_76-reaction_62;
	f_[7] = +reaction_44-reaction_48;
	f_[8] = +reaction_33-reaction_37;
	f_[9] = +reaction_1-reaction_2;
	f_[10] = +reaction_2;
	return GSL_SUCCESS;
}

/* ODE Jacobian: df(t,y;p)/dy   */
int AKAP79_jac(double t, const double y_[], double *jac_, double *dfdt_, void *par){
	double *p_=par;
	if (!y_ || !jac_) return 121;
/* 	constants   */
/* 	parameter values   */
	double kf_Rii_C__RiiP_C = p_[_kf_Rii_C__RiiP_C];                    /* [  0] */
	double kf_RiiP_CxcAMP__RiiP_C_cAMP = p_[_kf_RiiP_CxcAMP__RiiP_C_cAMP];       /* [  1] */
	double kf_RiiP_cAMPxC__RiiP_C_cAMP = p_[_kf_RiiP_cAMPxC__RiiP_C_cAMP];       /* [  2] */
	double kb_RiiP_cAMPxC__RiiP_C_cAMP = p_[_kb_RiiP_cAMPxC__RiiP_C_cAMP];       /* [  3] */
	double kb_RiiPXcAMP__RiiP_cAMP = p_[_kb_RiiPXcAMP__RiiP_cAMP];       /* [  4] */
	double kf_RiiPXcAMP__RiiP_cAMP = p_[_kf_RiiPXcAMP__RiiP_cAMP];       /* [  5] */
	double kf_RiiPxC__RiiP_C = p_[_kf_RiiPxC__RiiP_C];                  /* [  6] */
	double kb_RiiPxC__RiiP_C = p_[_kb_RiiPxC__RiiP_C];                  /* [  7] */
	double kf_cAMPxRii__Rii_cAMP = p_[_kf_cAMPxRii__Rii_cAMP];          /* [  8] */
	double kb_cAMPxRii__Rii_cAMP = p_[_kb_cAMPxRii__Rii_cAMP];          /* [  9] */
	double kf_Rii_CxcAMP__Rii_C_cAMP = p_[_kf_Rii_CxcAMP__Rii_C_cAMP];       /* [ 10] */
	double kb_Rii_CxcAMP__Rii_C_cAMP = p_[_kb_Rii_CxcAMP__Rii_C_cAMP];       /* [ 11] */
	double kf_RiixC__Rii_C = p_[_kf_RiixC__Rii_C];                      /* [ 12] */
	double kf_Rii_cAMPxC__Rii_C_cAMP = p_[_kf_Rii_cAMPxC__Rii_C_cAMP];       /* [ 13] */
	double kb_Rii_cAMPxC__Rii_C_cAMP = p_[_kb_Rii_cAMPxC__Rii_C_cAMP];       /* [ 14] */
	double kf_Rii_C_cAMP__RiiP_C_cAMP = p_[_kf_Rii_C_cAMP__RiiP_C_cAMP];       /* [ 15] */
	double kb_RiixC__Rii_C = p_[_kb_RiixC__Rii_C];                      /* [ 16] */
	double AKAPoff_1 = p_[_AKAPoff_1];                                  /* [ 17] */
	double AKAPoff_3 = p_[_AKAPoff_3];                                  /* [ 18] */
	double AKAPon_1 = p_[_AKAPon_1];                                    /* [ 19] */
	double AKAPon_3 = p_[_AKAPon_3];                                    /* [ 20] */
	double kf_C_AKAR4 = p_[_kf_C_AKAR4];                                /* [ 21] */
	double kb_C_AKAR4 = p_[_kb_C_AKAR4];                                /* [ 22] */
	double kcat_AKARp = p_[_kcat_AKARp];                                /* [ 23] */
	double kmOFF = p_[_kmOFF];                                          /* [ 24] */
	double kmON = p_[_kmON];                                            /* [ 25] */
	double KD_T = p_[_KD_T];                                            /* [ 26] */
	double b_AKAP = p_[_b_AKAP];                                        /* [ 27] */
	double AKAR4_ConservedConst = p_[_AKAR4_ConservedConst];            /* [ 28] */
	double CaN_ConservedConst = p_[_CaN_ConservedConst];                /* [ 29] */
	double Rii_C_ConservedConst = p_[_Rii_C_ConservedConst];            /* [ 30] */
	double cAMP_ConservedConst = p_[_cAMP_ConservedConst];              /* [ 31] */
	double Rii_ConservedConst = p_[_Rii_ConservedConst];                /* [ 32] */
/* 	state variables   */
	double RiiP = y_[_RiiP];                                            /* [  0] */
	double RiiP_cAMP = y_[_RiiP_cAMP];                                  /* [  1] */
	double RiiP_C = y_[_RiiP_C];                                        /* [  2] */
	double RiiP_C_cAMP = y_[_RiiP_C_cAMP];                              /* [  3] */
	double C = y_[_C];                                                  /* [  4] */
	double Rii_cAMP = y_[_Rii_cAMP];                                    /* [  5] */
	double Rii_C_cAMP = y_[_Rii_C_cAMP];                                /* [  6] */
	double RiiP_CaN = y_[_RiiP_CaN];                                    /* [  7] */
	double RiiP_cAMP_CaN = y_[_RiiP_cAMP_CaN];                          /* [  8] */
	double AKAR4_C = y_[_AKAR4_C];                                      /* [  9] */
	double AKAR4p = y_[_AKAR4p];                                        /* [ 10] */
/* 	expressions   */
	double kf_RiiP_cAMP_CaN__CaNXRii_cAMP = b_AKAP * AKAPon_1 + (1 - b_AKAP) * AKAPoff_1;
	double kb_RiiPxCaN__RiiP_CaN = b_AKAP*AKAPon_3  +  (1 - b_AKAP)* AKAPoff_3;
	double kf_RiiP_CaN__RiixCaN = b_AKAP * AKAPon_1 + (1 - b_AKAP) * AKAPoff_1;
	double kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN = b_AKAP*AKAPon_3  +  (1 - b_AKAP)* AKAPoff_3;
	double kf_RiiPxCaN__RiiP_CaN = b_AKAP * ((kb_RiiPxCaN__RiiP_CaN + kf_RiiP_cAMP_CaN__CaNXRii_cAMP)/kmON ) + (1 - b_AKAP) * (kb_RiiPxCaN__RiiP_CaN + kf_RiiP_cAMP_CaN__CaNXRii_cAMP) / kmOFF;
	double kf_CaNxRiiP_cAMP__RiiP_cAMP_CaN = b_AKAP * ((kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN + kf_RiiP_CaN__RiixCaN)/kmON) + (1 - b_AKAP) * (kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN + kf_RiiP_CaN__RiixCaN)/kmOFF;
	double kb_RiiP_CxcAMP__RiiP_C_cAMP = kf_RiiP_CxcAMP__RiiP_C_cAMP * KD_T;
	double AKAR4 = (AKAR4_ConservedConst - (AKAR4_C+AKAR4p));
	double CaN = (CaN_ConservedConst - (RiiP_CaN+RiiP_cAMP_CaN));
	double Rii_C = (Rii_C_ConservedConst - (RiiP_C+RiiP_C_cAMP+C+Rii_C_cAMP+AKAR4_C));
	double cAMP = (cAMP_ConservedConst - (RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN));
	double Rii = (Rii_ConservedConst - (RiiP+RiiP_cAMP+Rii_cAMP+RiiP_CaN+RiiP_cAMP_CaN-C-AKAR4_C));
	double reaction_51 = kf_Rii_C__RiiP_C*Rii_C;
	double reaction_14 = kf_RiiPxC__RiiP_C*RiiP*C - kb_RiiPxC__RiiP_C*RiiP_C;
	double reaction_12 = kf_RiiP_CxcAMP__RiiP_C_cAMP*RiiP_C*cAMP - kb_RiiP_CxcAMP__RiiP_C_cAMP*RiiP_C_cAMP;
	double reaction_43 = kf_RiiPXcAMP__RiiP_cAMP*cAMP*RiiP - kb_RiiPXcAMP__RiiP_cAMP*RiiP_cAMP;
	double reaction_23 = kf_RiiP_cAMPxC__RiiP_C_cAMP*RiiP_cAMP*C - kb_RiiP_cAMPxC__RiiP_C_cAMP*RiiP_C_cAMP;
	double reaction_78 = kf_cAMPxRii__Rii_cAMP*cAMP*Rii - kb_cAMPxRii__Rii_cAMP*Rii_cAMP;
	double reaction_56 = kf_Rii_CxcAMP__Rii_C_cAMP*Rii_C*cAMP - kb_Rii_CxcAMP__Rii_C_cAMP*Rii_C_cAMP;
	double reaction_76 = kf_Rii_cAMPxC__Rii_C_cAMP*Rii_cAMP*C - kb_Rii_cAMPxC__Rii_C_cAMP*Rii_C_cAMP;
	double reaction_62 = kf_Rii_C_cAMP__RiiP_C_cAMP*Rii_C_cAMP;
	double reaction_58 = kf_RiixC__Rii_C*Rii*C - kb_RiixC__Rii_C*Rii_C;
	double reaction_44 = kf_RiiPxCaN__RiiP_CaN*RiiP*CaN - kb_RiiPxCaN__RiiP_CaN*RiiP_CaN;
	double reaction_33 = kf_CaNxRiiP_cAMP__RiiP_cAMP_CaN*CaN*RiiP_cAMP - kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN*RiiP_cAMP_CaN;
	double reaction_48 = kf_RiiP_CaN__RiixCaN*RiiP_CaN;
	double reaction_37 = kf_RiiP_cAMP_CaN__CaNXRii_cAMP*RiiP_cAMP_CaN;
	double reaction_1 = kf_C_AKAR4*C*AKAR4 - kb_C_AKAR4*AKAR4_C;
	double reaction_2 = kcat_AKARp*AKAR4_C;
	memset(jac_,0,sizeof(double)*121); /* initialize with 0.0 */
	/*[ 0, 0]*/  jac_[0] = -(kf_RiiPxC__RiiP_C*C+kf_RiiPXcAMP__RiiP_cAMP*(cAMP_ConservedConst-(RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN))+((b_AKAP*(b_AKAP*AKAPon_3+(1-b_AKAP)*AKAPoff_3+b_AKAP*AKAPon_1+(1-b_AKAP)*AKAPoff_1))/kmON+((1-b_AKAP)*(b_AKAP*AKAPon_3+(1-b_AKAP)*AKAPoff_3+b_AKAP*AKAPon_1+(1-b_AKAP)*AKAPoff_1))/kmOFF)*(CaN_ConservedConst-(RiiP_CaN+RiiP_cAMP_CaN)));
	/*[ 0, 1]*/  jac_[1] = kf_RiiPXcAMP__RiiP_cAMP*RiiP+kb_RiiPXcAMP__RiiP_cAMP;
	/*[ 0, 2]*/  jac_[2] = kb_RiiPxC__RiiP_C;
	/*[ 0, 3]*/  jac_[3] = kf_RiiPXcAMP__RiiP_cAMP*RiiP;
	/*[ 0, 4]*/  jac_[4] = -kf_RiiPxC__RiiP_C*RiiP;
	/*[ 0, 5]*/  jac_[5] = kf_RiiPXcAMP__RiiP_cAMP*RiiP;
	/*[ 0, 6]*/  jac_[6] = kf_RiiPXcAMP__RiiP_cAMP*RiiP;
	/*[ 0, 7]*/  jac_[7] = ((b_AKAP*(b_AKAP*AKAPon_3+(1-b_AKAP)*AKAPoff_3+b_AKAP*AKAPon_1+(1-b_AKAP)*AKAPoff_1))/kmON+((1-b_AKAP)*(b_AKAP*AKAPon_3+(1-b_AKAP)*AKAPoff_3+b_AKAP*AKAPon_1+(1-b_AKAP)*AKAPoff_1))/kmOFF)*RiiP+b_AKAP*AKAPon_3+(1-b_AKAP)*AKAPoff_3;
	/*[ 0, 8]*/  jac_[8] = kf_RiiPXcAMP__RiiP_cAMP*RiiP+((b_AKAP*(b_AKAP*AKAPon_3+(1-b_AKAP)*AKAPoff_3+b_AKAP*AKAPon_1+(1-b_AKAP)*AKAPoff_1))/kmON+((1-b_AKAP)*(b_AKAP*AKAPon_3+(1-b_AKAP)*AKAPoff_3+b_AKAP*AKAPon_1+(1-b_AKAP)*AKAPoff_1))/kmOFF)*RiiP;
	/*[ 1, 0]*/  jac_[11] = kf_RiiPXcAMP__RiiP_cAMP*(cAMP_ConservedConst-(RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN));
	/*[ 1, 1]*/  jac_[12] = -(kf_RiiPXcAMP__RiiP_cAMP*RiiP+kb_RiiPXcAMP__RiiP_cAMP+kf_RiiP_cAMPxC__RiiP_C_cAMP*C+((b_AKAP*(b_AKAP*AKAPon_3+(1-b_AKAP)*AKAPoff_3+b_AKAP*AKAPon_1+(1-b_AKAP)*AKAPoff_1))/kmON+((1-b_AKAP)*(b_AKAP*AKAPon_3+(1-b_AKAP)*AKAPoff_3+b_AKAP*AKAPon_1+(1-b_AKAP)*AKAPoff_1))/kmOFF)*(CaN_ConservedConst-(RiiP_CaN+RiiP_cAMP_CaN)));
	/*[ 1, 3]*/  jac_[14] = kb_RiiP_cAMPxC__RiiP_C_cAMP-kf_RiiPXcAMP__RiiP_cAMP*RiiP;
	/*[ 1, 4]*/  jac_[15] = -kf_RiiP_cAMPxC__RiiP_C_cAMP*RiiP_cAMP;
	/*[ 1, 5]*/  jac_[16] = -kf_RiiPXcAMP__RiiP_cAMP*RiiP;
	/*[ 1, 6]*/  jac_[17] = -kf_RiiPXcAMP__RiiP_cAMP*RiiP;
	/*[ 1, 7]*/  jac_[18] = ((b_AKAP*(b_AKAP*AKAPon_3+(1-b_AKAP)*AKAPoff_3+b_AKAP*AKAPon_1+(1-b_AKAP)*AKAPoff_1))/kmON+((1-b_AKAP)*(b_AKAP*AKAPon_3+(1-b_AKAP)*AKAPoff_3+b_AKAP*AKAPon_1+(1-b_AKAP)*AKAPoff_1))/kmOFF)*RiiP_cAMP;
	/*[ 1, 8]*/  jac_[19] = ((b_AKAP*(b_AKAP*AKAPon_3+(1-b_AKAP)*AKAPoff_3+b_AKAP*AKAPon_1+(1-b_AKAP)*AKAPoff_1))/kmON+((1-b_AKAP)*(b_AKAP*AKAPon_3+(1-b_AKAP)*AKAPoff_3+b_AKAP*AKAPon_1+(1-b_AKAP)*AKAPoff_1))/kmOFF)*RiiP_cAMP+b_AKAP*AKAPon_3+(1-b_AKAP)*AKAPoff_3-kf_RiiPXcAMP__RiiP_cAMP*RiiP;
	/*[ 2, 0]*/  jac_[22] = kf_RiiPxC__RiiP_C*C;
	/*[ 2, 1]*/  jac_[23] = kf_RiiP_CxcAMP__RiiP_C_cAMP*RiiP_C;
	/*[ 2, 2]*/  jac_[24] = -(kb_RiiPxC__RiiP_C+kf_Rii_C__RiiP_C+kf_RiiP_CxcAMP__RiiP_C_cAMP*(cAMP_ConservedConst-(RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN)));
	/*[ 2, 3]*/  jac_[25] = kf_RiiP_CxcAMP__RiiP_C_cAMP*RiiP_C+kf_RiiP_CxcAMP__RiiP_C_cAMP*KD_T-kf_Rii_C__RiiP_C;
	/*[ 2, 4]*/  jac_[26] = kf_RiiPxC__RiiP_C*RiiP-kf_Rii_C__RiiP_C;
	/*[ 2, 5]*/  jac_[27] = kf_RiiP_CxcAMP__RiiP_C_cAMP*RiiP_C;
	/*[ 2, 6]*/  jac_[28] = kf_RiiP_CxcAMP__RiiP_C_cAMP*RiiP_C-kf_Rii_C__RiiP_C;
	/*[ 2, 8]*/  jac_[30] = kf_RiiP_CxcAMP__RiiP_C_cAMP*RiiP_C;
	/*[ 2, 9]*/  jac_[31] = -kf_Rii_C__RiiP_C;
	/*[ 3, 1]*/  jac_[34] = kf_RiiP_cAMPxC__RiiP_C_cAMP*C-kf_RiiP_CxcAMP__RiiP_C_cAMP*RiiP_C;
	/*[ 3, 2]*/  jac_[35] = kf_RiiP_CxcAMP__RiiP_C_cAMP*(cAMP_ConservedConst-(RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN));
	/*[ 3, 3]*/  jac_[36] = -(kb_RiiP_cAMPxC__RiiP_C_cAMP+kf_RiiP_CxcAMP__RiiP_C_cAMP*RiiP_C+kf_RiiP_CxcAMP__RiiP_C_cAMP*KD_T);
	/*[ 3, 4]*/  jac_[37] = kf_RiiP_cAMPxC__RiiP_C_cAMP*RiiP_cAMP;
	/*[ 3, 5]*/  jac_[38] = -kf_RiiP_CxcAMP__RiiP_C_cAMP*RiiP_C;
	/*[ 3, 6]*/  jac_[39] = kf_Rii_C_cAMP__RiiP_C_cAMP-kf_RiiP_CxcAMP__RiiP_C_cAMP*RiiP_C;
	/*[ 3, 8]*/  jac_[41] = -kf_RiiP_CxcAMP__RiiP_C_cAMP*RiiP_C;
	/*[ 4, 0]*/  jac_[44] = kf_RiixC__Rii_C*C-kf_RiiPxC__RiiP_C*C;
	/*[ 4, 1]*/  jac_[45] = kf_RiixC__Rii_C*C-kf_RiiP_cAMPxC__RiiP_C_cAMP*C;
	/*[ 4, 2]*/  jac_[46] = kb_RiiPxC__RiiP_C-kb_RiixC__Rii_C;
	/*[ 4, 3]*/  jac_[47] = kb_RiiP_cAMPxC__RiiP_C_cAMP-kb_RiixC__Rii_C;
	/*[ 4, 4]*/  jac_[48] = -(kf_RiiPxC__RiiP_C*RiiP+kf_RiiP_cAMPxC__RiiP_C_cAMP*RiiP_cAMP+kf_Rii_cAMPxC__Rii_C_cAMP*Rii_cAMP+kf_RiixC__Rii_C*(Rii_ConservedConst-(RiiP+RiiP_cAMP+Rii_cAMP+RiiP_CaN+RiiP_cAMP_CaN-C-AKAR4_C))+kf_RiixC__Rii_C*C+kb_RiixC__Rii_C+kf_C_AKAR4*(AKAR4_ConservedConst-(AKAR4_C+AKAR4p)));
	/*[ 4, 5]*/  jac_[49] = kf_RiixC__Rii_C*C-kf_Rii_cAMPxC__Rii_C_cAMP*C;
	/*[ 4, 6]*/  jac_[50] = kb_Rii_cAMPxC__Rii_C_cAMP-kb_RiixC__Rii_C;
	/*[ 4, 7]*/  jac_[51] = kf_RiixC__Rii_C*C;
	/*[ 4, 8]*/  jac_[52] = kf_RiixC__Rii_C*C;
	/*[ 4, 9]*/  jac_[53] = kf_C_AKAR4*C+kb_C_AKAR4-(kf_RiixC__Rii_C*C+kb_RiixC__Rii_C)+kcat_AKARp;
	/*[ 4,10]*/  jac_[54] = kf_C_AKAR4*C;
	/*[ 5, 0]*/  jac_[55] = -kf_cAMPxRii__Rii_cAMP*(cAMP_ConservedConst-(RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN));
	/*[ 5, 1]*/  jac_[56] = -(kf_cAMPxRii__Rii_cAMP*(Rii_ConservedConst-(RiiP+RiiP_cAMP+Rii_cAMP+RiiP_CaN+RiiP_cAMP_CaN-C-AKAR4_C))+kf_cAMPxRii__Rii_cAMP*(cAMP_ConservedConst-(RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN)));
	/*[ 5, 3]*/  jac_[58] = -kf_cAMPxRii__Rii_cAMP*(Rii_ConservedConst-(RiiP+RiiP_cAMP+Rii_cAMP+RiiP_CaN+RiiP_cAMP_CaN-C-AKAR4_C));
	/*[ 5, 4]*/  jac_[59] = kf_cAMPxRii__Rii_cAMP*(cAMP_ConservedConst-(RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN))-kf_Rii_cAMPxC__Rii_C_cAMP*Rii_cAMP;
	/*[ 5, 5]*/  jac_[60] = -(kf_cAMPxRii__Rii_cAMP*(Rii_ConservedConst-(RiiP+RiiP_cAMP+Rii_cAMP+RiiP_CaN+RiiP_cAMP_CaN-C-AKAR4_C))+kf_cAMPxRii__Rii_cAMP*(cAMP_ConservedConst-(RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN))+kb_cAMPxRii__Rii_cAMP+kf_Rii_cAMPxC__Rii_C_cAMP*C);
	/*[ 5, 6]*/  jac_[61] = kb_Rii_cAMPxC__Rii_C_cAMP-kf_cAMPxRii__Rii_cAMP*(Rii_ConservedConst-(RiiP+RiiP_cAMP+Rii_cAMP+RiiP_CaN+RiiP_cAMP_CaN-C-AKAR4_C));
	/*[ 5, 7]*/  jac_[62] = -kf_cAMPxRii__Rii_cAMP*(cAMP_ConservedConst-(RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN));
	/*[ 5, 8]*/  jac_[63] = b_AKAP*AKAPon_1+(1-b_AKAP)*AKAPoff_1-(kf_cAMPxRii__Rii_cAMP*(Rii_ConservedConst-(RiiP+RiiP_cAMP+Rii_cAMP+RiiP_CaN+RiiP_cAMP_CaN-C-AKAR4_C))+kf_cAMPxRii__Rii_cAMP*(cAMP_ConservedConst-(RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN)));
	/*[ 5, 9]*/  jac_[64] = kf_cAMPxRii__Rii_cAMP*(cAMP_ConservedConst-(RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN));
	/*[ 6, 1]*/  jac_[67] = -kf_Rii_CxcAMP__Rii_C_cAMP*(Rii_C_ConservedConst-(RiiP_C+RiiP_C_cAMP+C+Rii_C_cAMP+AKAR4_C));
	/*[ 6, 2]*/  jac_[68] = -kf_Rii_CxcAMP__Rii_C_cAMP*(cAMP_ConservedConst-(RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN));
	/*[ 6, 3]*/  jac_[69] = -(kf_Rii_CxcAMP__Rii_C_cAMP*(cAMP_ConservedConst-(RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN))+kf_Rii_CxcAMP__Rii_C_cAMP*(Rii_C_ConservedConst-(RiiP_C+RiiP_C_cAMP+C+Rii_C_cAMP+AKAR4_C)));
	/*[ 6, 4]*/  jac_[70] = kf_Rii_cAMPxC__Rii_C_cAMP*Rii_cAMP-kf_Rii_CxcAMP__Rii_C_cAMP*(cAMP_ConservedConst-(RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN));
	/*[ 6, 5]*/  jac_[71] = kf_Rii_cAMPxC__Rii_C_cAMP*C-kf_Rii_CxcAMP__Rii_C_cAMP*(Rii_C_ConservedConst-(RiiP_C+RiiP_C_cAMP+C+Rii_C_cAMP+AKAR4_C));
	/*[ 6, 6]*/  jac_[72] = -(kb_Rii_cAMPxC__Rii_C_cAMP+kf_Rii_CxcAMP__Rii_C_cAMP*(cAMP_ConservedConst-(RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN))+kf_Rii_CxcAMP__Rii_C_cAMP*(Rii_C_ConservedConst-(RiiP_C+RiiP_C_cAMP+C+Rii_C_cAMP+AKAR4_C))+kb_Rii_CxcAMP__Rii_C_cAMP+kf_Rii_C_cAMP__RiiP_C_cAMP);
	/*[ 6, 8]*/  jac_[74] = -kf_Rii_CxcAMP__Rii_C_cAMP*(Rii_C_ConservedConst-(RiiP_C+RiiP_C_cAMP+C+Rii_C_cAMP+AKAR4_C));
	/*[ 6, 9]*/  jac_[75] = -kf_Rii_CxcAMP__Rii_C_cAMP*(cAMP_ConservedConst-(RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN));
	/*[ 7, 0]*/  jac_[77] = ((b_AKAP*(b_AKAP*AKAPon_3+(1-b_AKAP)*AKAPoff_3+b_AKAP*AKAPon_1+(1-b_AKAP)*AKAPoff_1))/kmON+((1-b_AKAP)*(b_AKAP*AKAPon_3+(1-b_AKAP)*AKAPoff_3+b_AKAP*AKAPon_1+(1-b_AKAP)*AKAPoff_1))/kmOFF)*(CaN_ConservedConst-(RiiP_CaN+RiiP_cAMP_CaN));
	/*[ 7, 7]*/  jac_[84] = -(((b_AKAP*(b_AKAP*AKAPon_3+(1-b_AKAP)*AKAPoff_3+b_AKAP*AKAPon_1+(1-b_AKAP)*AKAPoff_1))/kmON+((1-b_AKAP)*(b_AKAP*AKAPon_3+(1-b_AKAP)*AKAPoff_3+b_AKAP*AKAPon_1+(1-b_AKAP)*AKAPoff_1))/kmOFF)*RiiP+b_AKAP*AKAPon_3+(1-b_AKAP)*AKAPoff_3+b_AKAP*AKAPon_1+(1-b_AKAP)*AKAPoff_1);
	/*[ 7, 8]*/  jac_[85] = -((b_AKAP*(b_AKAP*AKAPon_3+(1-b_AKAP)*AKAPoff_3+b_AKAP*AKAPon_1+(1-b_AKAP)*AKAPoff_1))/kmON+((1-b_AKAP)*(b_AKAP*AKAPon_3+(1-b_AKAP)*AKAPoff_3+b_AKAP*AKAPon_1+(1-b_AKAP)*AKAPoff_1))/kmOFF)*RiiP;
	/*[ 8, 1]*/  jac_[89] = ((b_AKAP*(b_AKAP*AKAPon_3+(1-b_AKAP)*AKAPoff_3+b_AKAP*AKAPon_1+(1-b_AKAP)*AKAPoff_1))/kmON+((1-b_AKAP)*(b_AKAP*AKAPon_3+(1-b_AKAP)*AKAPoff_3+b_AKAP*AKAPon_1+(1-b_AKAP)*AKAPoff_1))/kmOFF)*(CaN_ConservedConst-(RiiP_CaN+RiiP_cAMP_CaN));
	/*[ 8, 7]*/  jac_[95] = -((b_AKAP*(b_AKAP*AKAPon_3+(1-b_AKAP)*AKAPoff_3+b_AKAP*AKAPon_1+(1-b_AKAP)*AKAPoff_1))/kmON+((1-b_AKAP)*(b_AKAP*AKAPon_3+(1-b_AKAP)*AKAPoff_3+b_AKAP*AKAPon_1+(1-b_AKAP)*AKAPoff_1))/kmOFF)*RiiP_cAMP;
	/*[ 8, 8]*/  jac_[96] = -(((b_AKAP*(b_AKAP*AKAPon_3+(1-b_AKAP)*AKAPoff_3+b_AKAP*AKAPon_1+(1-b_AKAP)*AKAPoff_1))/kmON+((1-b_AKAP)*(b_AKAP*AKAPon_3+(1-b_AKAP)*AKAPoff_3+b_AKAP*AKAPon_1+(1-b_AKAP)*AKAPoff_1))/kmOFF)*RiiP_cAMP+b_AKAP*AKAPon_3+(1-b_AKAP)*AKAPoff_3+b_AKAP*AKAPon_1+(1-b_AKAP)*AKAPoff_1);
	/*[ 9, 4]*/  jac_[103] = kf_C_AKAR4*(AKAR4_ConservedConst-(AKAR4_C+AKAR4p));
	/*[ 9, 9]*/  jac_[108] = -(kf_C_AKAR4*C+kb_C_AKAR4+kcat_AKARp);
	/*[ 9,10]*/  jac_[109] = -kf_C_AKAR4*C;
	/*[10, 9]*/  jac_[119] = kcat_AKARp;
	return GSL_SUCCESS;
}

/* ODE parameter Jacobian: df(t,y;p)/dp   */
int AKAP79_jacp(double t, const double y_[], double *jacp_, double *dfdt_, void *par){
	double *p_=par;
	if (!y_ || !jacp_) return 363;
/* 	constants   */
/* 	parameter values   */
	double kf_Rii_C__RiiP_C = p_[_kf_Rii_C__RiiP_C];                    /* [  0] */
	double kf_RiiP_CxcAMP__RiiP_C_cAMP = p_[_kf_RiiP_CxcAMP__RiiP_C_cAMP];       /* [  1] */
	double kf_RiiP_cAMPxC__RiiP_C_cAMP = p_[_kf_RiiP_cAMPxC__RiiP_C_cAMP];       /* [  2] */
	double kb_RiiP_cAMPxC__RiiP_C_cAMP = p_[_kb_RiiP_cAMPxC__RiiP_C_cAMP];       /* [  3] */
	double kb_RiiPXcAMP__RiiP_cAMP = p_[_kb_RiiPXcAMP__RiiP_cAMP];       /* [  4] */
	double kf_RiiPXcAMP__RiiP_cAMP = p_[_kf_RiiPXcAMP__RiiP_cAMP];       /* [  5] */
	double kf_RiiPxC__RiiP_C = p_[_kf_RiiPxC__RiiP_C];                  /* [  6] */
	double kb_RiiPxC__RiiP_C = p_[_kb_RiiPxC__RiiP_C];                  /* [  7] */
	double kf_cAMPxRii__Rii_cAMP = p_[_kf_cAMPxRii__Rii_cAMP];          /* [  8] */
	double kb_cAMPxRii__Rii_cAMP = p_[_kb_cAMPxRii__Rii_cAMP];          /* [  9] */
	double kf_Rii_CxcAMP__Rii_C_cAMP = p_[_kf_Rii_CxcAMP__Rii_C_cAMP];       /* [ 10] */
	double kb_Rii_CxcAMP__Rii_C_cAMP = p_[_kb_Rii_CxcAMP__Rii_C_cAMP];       /* [ 11] */
	double kf_RiixC__Rii_C = p_[_kf_RiixC__Rii_C];                      /* [ 12] */
	double kf_Rii_cAMPxC__Rii_C_cAMP = p_[_kf_Rii_cAMPxC__Rii_C_cAMP];       /* [ 13] */
	double kb_Rii_cAMPxC__Rii_C_cAMP = p_[_kb_Rii_cAMPxC__Rii_C_cAMP];       /* [ 14] */
	double kf_Rii_C_cAMP__RiiP_C_cAMP = p_[_kf_Rii_C_cAMP__RiiP_C_cAMP];       /* [ 15] */
	double kb_RiixC__Rii_C = p_[_kb_RiixC__Rii_C];                      /* [ 16] */
	double AKAPoff_1 = p_[_AKAPoff_1];                                  /* [ 17] */
	double AKAPoff_3 = p_[_AKAPoff_3];                                  /* [ 18] */
	double AKAPon_1 = p_[_AKAPon_1];                                    /* [ 19] */
	double AKAPon_3 = p_[_AKAPon_3];                                    /* [ 20] */
	double kf_C_AKAR4 = p_[_kf_C_AKAR4];                                /* [ 21] */
	double kb_C_AKAR4 = p_[_kb_C_AKAR4];                                /* [ 22] */
	double kcat_AKARp = p_[_kcat_AKARp];                                /* [ 23] */
	double kmOFF = p_[_kmOFF];                                          /* [ 24] */
	double kmON = p_[_kmON];                                            /* [ 25] */
	double KD_T = p_[_KD_T];                                            /* [ 26] */
	double b_AKAP = p_[_b_AKAP];                                        /* [ 27] */
	double AKAR4_ConservedConst = p_[_AKAR4_ConservedConst];            /* [ 28] */
	double CaN_ConservedConst = p_[_CaN_ConservedConst];                /* [ 29] */
	double Rii_C_ConservedConst = p_[_Rii_C_ConservedConst];            /* [ 30] */
	double cAMP_ConservedConst = p_[_cAMP_ConservedConst];              /* [ 31] */
	double Rii_ConservedConst = p_[_Rii_ConservedConst];                /* [ 32] */
/* 	state variables   */
	double RiiP = y_[_RiiP];                                            /* [  0] */
	double RiiP_cAMP = y_[_RiiP_cAMP];                                  /* [  1] */
	double RiiP_C = y_[_RiiP_C];                                        /* [  2] */
	double RiiP_C_cAMP = y_[_RiiP_C_cAMP];                              /* [  3] */
	double C = y_[_C];                                                  /* [  4] */
	double Rii_cAMP = y_[_Rii_cAMP];                                    /* [  5] */
	double Rii_C_cAMP = y_[_Rii_C_cAMP];                                /* [  6] */
	double RiiP_CaN = y_[_RiiP_CaN];                                    /* [  7] */
	double RiiP_cAMP_CaN = y_[_RiiP_cAMP_CaN];                          /* [  8] */
	double AKAR4_C = y_[_AKAR4_C];                                      /* [  9] */
	double AKAR4p = y_[_AKAR4p];                                        /* [ 10] */
/* 	expressions   */
	double kf_RiiP_cAMP_CaN__CaNXRii_cAMP = b_AKAP * AKAPon_1 + (1 - b_AKAP) * AKAPoff_1;
	double kb_RiiPxCaN__RiiP_CaN = b_AKAP*AKAPon_3  +  (1 - b_AKAP)* AKAPoff_3;
	double kf_RiiP_CaN__RiixCaN = b_AKAP * AKAPon_1 + (1 - b_AKAP) * AKAPoff_1;
	double kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN = b_AKAP*AKAPon_3  +  (1 - b_AKAP)* AKAPoff_3;
	double kf_RiiPxCaN__RiiP_CaN = b_AKAP * ((kb_RiiPxCaN__RiiP_CaN + kf_RiiP_cAMP_CaN__CaNXRii_cAMP)/kmON ) + (1 - b_AKAP) * (kb_RiiPxCaN__RiiP_CaN + kf_RiiP_cAMP_CaN__CaNXRii_cAMP) / kmOFF;
	double kf_CaNxRiiP_cAMP__RiiP_cAMP_CaN = b_AKAP * ((kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN + kf_RiiP_CaN__RiixCaN)/kmON) + (1 - b_AKAP) * (kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN + kf_RiiP_CaN__RiixCaN)/kmOFF;
	double kb_RiiP_CxcAMP__RiiP_C_cAMP = kf_RiiP_CxcAMP__RiiP_C_cAMP * KD_T;
	double AKAR4 = (AKAR4_ConservedConst - (AKAR4_C+AKAR4p));
	double CaN = (CaN_ConservedConst - (RiiP_CaN+RiiP_cAMP_CaN));
	double Rii_C = (Rii_C_ConservedConst - (RiiP_C+RiiP_C_cAMP+C+Rii_C_cAMP+AKAR4_C));
	double cAMP = (cAMP_ConservedConst - (RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN));
	double Rii = (Rii_ConservedConst - (RiiP+RiiP_cAMP+Rii_cAMP+RiiP_CaN+RiiP_cAMP_CaN-C-AKAR4_C));
	double reaction_51 = kf_Rii_C__RiiP_C*Rii_C;
	double reaction_14 = kf_RiiPxC__RiiP_C*RiiP*C - kb_RiiPxC__RiiP_C*RiiP_C;
	double reaction_12 = kf_RiiP_CxcAMP__RiiP_C_cAMP*RiiP_C*cAMP - kb_RiiP_CxcAMP__RiiP_C_cAMP*RiiP_C_cAMP;
	double reaction_43 = kf_RiiPXcAMP__RiiP_cAMP*cAMP*RiiP - kb_RiiPXcAMP__RiiP_cAMP*RiiP_cAMP;
	double reaction_23 = kf_RiiP_cAMPxC__RiiP_C_cAMP*RiiP_cAMP*C - kb_RiiP_cAMPxC__RiiP_C_cAMP*RiiP_C_cAMP;
	double reaction_78 = kf_cAMPxRii__Rii_cAMP*cAMP*Rii - kb_cAMPxRii__Rii_cAMP*Rii_cAMP;
	double reaction_56 = kf_Rii_CxcAMP__Rii_C_cAMP*Rii_C*cAMP - kb_Rii_CxcAMP__Rii_C_cAMP*Rii_C_cAMP;
	double reaction_76 = kf_Rii_cAMPxC__Rii_C_cAMP*Rii_cAMP*C - kb_Rii_cAMPxC__Rii_C_cAMP*Rii_C_cAMP;
	double reaction_62 = kf_Rii_C_cAMP__RiiP_C_cAMP*Rii_C_cAMP;
	double reaction_58 = kf_RiixC__Rii_C*Rii*C - kb_RiixC__Rii_C*Rii_C;
	double reaction_44 = kf_RiiPxCaN__RiiP_CaN*RiiP*CaN - kb_RiiPxCaN__RiiP_CaN*RiiP_CaN;
	double reaction_33 = kf_CaNxRiiP_cAMP__RiiP_cAMP_CaN*CaN*RiiP_cAMP - kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN*RiiP_cAMP_CaN;
	double reaction_48 = kf_RiiP_CaN__RiixCaN*RiiP_CaN;
	double reaction_37 = kf_RiiP_cAMP_CaN__CaNXRii_cAMP*RiiP_cAMP_CaN;
	double reaction_1 = kf_C_AKAR4*C*AKAR4 - kb_C_AKAR4*AKAR4_C;
	double reaction_2 = kcat_AKARp*AKAR4_C;
	memset(jacp_,0,sizeof(double)*363); /* initialize with 0.0 */
	/*[ 0, 4]*/  jacp_[4] = RiiP_cAMP;
	/*[ 0, 5]*/  jacp_[5] = -(cAMP_ConservedConst-(RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN))*RiiP;
	/*[ 0, 6]*/  jacp_[6] = -RiiP*C;
	/*[ 0, 7]*/  jacp_[7] = RiiP_C;
	/*[ 0,17]*/  jacp_[17] = -((b_AKAP*(1-b_AKAP))/kmON+gsl_pow_2(1-b_AKAP)/kmOFF)*RiiP*(CaN_ConservedConst-(RiiP_CaN+RiiP_cAMP_CaN));
	/*[ 0,18]*/  jacp_[18] = (1-b_AKAP)*RiiP_CaN-((b_AKAP*(1-b_AKAP))/kmON+gsl_pow_2(1-b_AKAP)/kmOFF)*RiiP*(CaN_ConservedConst-(RiiP_CaN+RiiP_cAMP_CaN));
	/*[ 0,19]*/  jacp_[19] = -(gsl_pow_2(b_AKAP)/kmON+((1-b_AKAP)*b_AKAP)/kmOFF)*RiiP*(CaN_ConservedConst-(RiiP_CaN+RiiP_cAMP_CaN));
	/*[ 0,20]*/  jacp_[20] = b_AKAP*RiiP_CaN-(gsl_pow_2(b_AKAP)/kmON+((1-b_AKAP)*b_AKAP)/kmOFF)*RiiP*(CaN_ConservedConst-(RiiP_CaN+RiiP_cAMP_CaN));
	/*[ 0,24]*/  jacp_[24] = ((CaN_ConservedConst-(RiiP_CaN+RiiP_cAMP_CaN))*RiiP*(1-b_AKAP)*(b_AKAP*AKAPon_3+(1-b_AKAP)*AKAPoff_3+b_AKAP*AKAPon_1+(1-b_AKAP)*AKAPoff_1))/gsl_pow_2(kmOFF);
	/*[ 0,25]*/  jacp_[25] = ((CaN_ConservedConst-(RiiP_CaN+RiiP_cAMP_CaN))*RiiP*b_AKAP*(b_AKAP*AKAPon_3+(1-b_AKAP)*AKAPoff_3+b_AKAP*AKAPon_1+(1-b_AKAP)*AKAPoff_1))/gsl_pow_2(kmON);
	/*[ 0,27]*/  jacp_[27] = (AKAPon_3-AKAPoff_3)*RiiP_CaN-((b_AKAP*(AKAPon_3-AKAPoff_3+AKAPon_1-AKAPoff_1)+b_AKAP*AKAPon_3+(1-b_AKAP)*AKAPoff_3+b_AKAP*AKAPon_1+(1-b_AKAP)*AKAPoff_1)/kmON+((1-b_AKAP)*(AKAPon_3-AKAPoff_3+AKAPon_1-AKAPoff_1)-(b_AKAP*AKAPon_3+(1-b_AKAP)*AKAPoff_3+b_AKAP*AKAPon_1+(1-b_AKAP)*AKAPoff_1))/kmOFF)*RiiP*(CaN_ConservedConst-(RiiP_CaN+RiiP_cAMP_CaN));
	/*[ 0,29]*/  jacp_[29] = -((b_AKAP*(b_AKAP*AKAPon_3+(1-b_AKAP)*AKAPoff_3+b_AKAP*AKAPon_1+(1-b_AKAP)*AKAPoff_1))/kmON+((1-b_AKAP)*(b_AKAP*AKAPon_3+(1-b_AKAP)*AKAPoff_3+b_AKAP*AKAPon_1+(1-b_AKAP)*AKAPoff_1))/kmOFF)*RiiP;
	/*[ 0,31]*/  jacp_[31] = -kf_RiiPXcAMP__RiiP_cAMP*RiiP;
	/*[ 1, 2]*/  jacp_[35] = -RiiP_cAMP*C;
	/*[ 1, 3]*/  jacp_[36] = RiiP_C_cAMP;
	/*[ 1, 4]*/  jacp_[37] = -RiiP_cAMP;
	/*[ 1, 5]*/  jacp_[38] = (cAMP_ConservedConst-(RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN))*RiiP;
	/*[ 1,17]*/  jacp_[50] = -((b_AKAP*(1-b_AKAP))/kmON+gsl_pow_2(1-b_AKAP)/kmOFF)*(CaN_ConservedConst-(RiiP_CaN+RiiP_cAMP_CaN))*RiiP_cAMP;
	/*[ 1,18]*/  jacp_[51] = (1-b_AKAP)*RiiP_cAMP_CaN-((b_AKAP*(1-b_AKAP))/kmON+gsl_pow_2(1-b_AKAP)/kmOFF)*(CaN_ConservedConst-(RiiP_CaN+RiiP_cAMP_CaN))*RiiP_cAMP;
	/*[ 1,19]*/  jacp_[52] = -(gsl_pow_2(b_AKAP)/kmON+((1-b_AKAP)*b_AKAP)/kmOFF)*(CaN_ConservedConst-(RiiP_CaN+RiiP_cAMP_CaN))*RiiP_cAMP;
	/*[ 1,20]*/  jacp_[53] = b_AKAP*RiiP_cAMP_CaN-(gsl_pow_2(b_AKAP)/kmON+((1-b_AKAP)*b_AKAP)/kmOFF)*(CaN_ConservedConst-(RiiP_CaN+RiiP_cAMP_CaN))*RiiP_cAMP;
	/*[ 1,24]*/  jacp_[57] = (RiiP_cAMP*(CaN_ConservedConst-(RiiP_CaN+RiiP_cAMP_CaN))*(1-b_AKAP)*(b_AKAP*AKAPon_3+(1-b_AKAP)*AKAPoff_3+b_AKAP*AKAPon_1+(1-b_AKAP)*AKAPoff_1))/gsl_pow_2(kmOFF);
	/*[ 1,25]*/  jacp_[58] = (RiiP_cAMP*(CaN_ConservedConst-(RiiP_CaN+RiiP_cAMP_CaN))*b_AKAP*(b_AKAP*AKAPon_3+(1-b_AKAP)*AKAPoff_3+b_AKAP*AKAPon_1+(1-b_AKAP)*AKAPoff_1))/gsl_pow_2(kmON);
	/*[ 1,27]*/  jacp_[60] = (AKAPon_3-AKAPoff_3)*RiiP_cAMP_CaN-((b_AKAP*(AKAPon_3-AKAPoff_3+AKAPon_1-AKAPoff_1)+b_AKAP*AKAPon_3+(1-b_AKAP)*AKAPoff_3+b_AKAP*AKAPon_1+(1-b_AKAP)*AKAPoff_1)/kmON+((1-b_AKAP)*(AKAPon_3-AKAPoff_3+AKAPon_1-AKAPoff_1)-(b_AKAP*AKAPon_3+(1-b_AKAP)*AKAPoff_3+b_AKAP*AKAPon_1+(1-b_AKAP)*AKAPoff_1))/kmOFF)*(CaN_ConservedConst-(RiiP_CaN+RiiP_cAMP_CaN))*RiiP_cAMP;
	/*[ 1,29]*/  jacp_[62] = -((b_AKAP*(b_AKAP*AKAPon_3+(1-b_AKAP)*AKAPoff_3+b_AKAP*AKAPon_1+(1-b_AKAP)*AKAPoff_1))/kmON+((1-b_AKAP)*(b_AKAP*AKAPon_3+(1-b_AKAP)*AKAPoff_3+b_AKAP*AKAPon_1+(1-b_AKAP)*AKAPoff_1))/kmOFF)*RiiP_cAMP;
	/*[ 1,31]*/  jacp_[64] = kf_RiiPXcAMP__RiiP_cAMP*RiiP;
	/*[ 2, 0]*/  jacp_[66] = Rii_C_ConservedConst-(RiiP_C+RiiP_C_cAMP+C+Rii_C_cAMP+AKAR4_C);
	/*[ 2, 1]*/  jacp_[67] = KD_T*RiiP_C_cAMP-RiiP_C*(cAMP_ConservedConst-(RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN));
	/*[ 2, 6]*/  jacp_[72] = RiiP*C;
	/*[ 2, 7]*/  jacp_[73] = -RiiP_C;
	/*[ 2,26]*/  jacp_[92] = kf_RiiP_CxcAMP__RiiP_C_cAMP*RiiP_C_cAMP;
	/*[ 2,30]*/  jacp_[96] = kf_Rii_C__RiiP_C;
	/*[ 2,31]*/  jacp_[97] = -kf_RiiP_CxcAMP__RiiP_C_cAMP*RiiP_C;
	/*[ 3, 1]*/  jacp_[100] = RiiP_C*(cAMP_ConservedConst-(RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN))-KD_T*RiiP_C_cAMP;
	/*[ 3, 2]*/  jacp_[101] = RiiP_cAMP*C;
	/*[ 3, 3]*/  jacp_[102] = -RiiP_C_cAMP;
	/*[ 3,15]*/  jacp_[114] = Rii_C_cAMP;
	/*[ 3,26]*/  jacp_[125] = -kf_RiiP_CxcAMP__RiiP_C_cAMP*RiiP_C_cAMP;
	/*[ 3,31]*/  jacp_[130] = kf_RiiP_CxcAMP__RiiP_C_cAMP*RiiP_C;
	/*[ 4, 2]*/  jacp_[134] = -RiiP_cAMP*C;
	/*[ 4, 3]*/  jacp_[135] = RiiP_C_cAMP;
	/*[ 4, 6]*/  jacp_[138] = -RiiP*C;
	/*[ 4, 7]*/  jacp_[139] = RiiP_C;
	/*[ 4,12]*/  jacp_[144] = -(Rii_ConservedConst-(RiiP+RiiP_cAMP+Rii_cAMP+RiiP_CaN+RiiP_cAMP_CaN-C-AKAR4_C))*C;
	/*[ 4,13]*/  jacp_[145] = -Rii_cAMP*C;
	/*[ 4,14]*/  jacp_[146] = Rii_C_cAMP;
	/*[ 4,16]*/  jacp_[148] = Rii_C_ConservedConst-(RiiP_C+RiiP_C_cAMP+C+Rii_C_cAMP+AKAR4_C);
	/*[ 4,21]*/  jacp_[153] = -C*(AKAR4_ConservedConst-(AKAR4_C+AKAR4p));
	/*[ 4,22]*/  jacp_[154] = AKAR4_C;
	/*[ 4,23]*/  jacp_[155] = AKAR4_C;
	/*[ 4,28]*/  jacp_[160] = -kf_C_AKAR4*C;
	/*[ 4,30]*/  jacp_[162] = kb_RiixC__Rii_C;
	/*[ 4,32]*/  jacp_[164] = -kf_RiixC__Rii_C*C;
	/*[ 5, 8]*/  jacp_[173] = (cAMP_ConservedConst-(RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN))*(Rii_ConservedConst-(RiiP+RiiP_cAMP+Rii_cAMP+RiiP_CaN+RiiP_cAMP_CaN-C-AKAR4_C));
	/*[ 5, 9]*/  jacp_[174] = -Rii_cAMP;
	/*[ 5,13]*/  jacp_[178] = -Rii_cAMP*C;
	/*[ 5,14]*/  jacp_[179] = Rii_C_cAMP;
	/*[ 5,17]*/  jacp_[182] = (1-b_AKAP)*RiiP_cAMP_CaN;
	/*[ 5,19]*/  jacp_[184] = b_AKAP*RiiP_cAMP_CaN;
	/*[ 5,27]*/  jacp_[192] = (AKAPon_1-AKAPoff_1)*RiiP_cAMP_CaN;
	/*[ 5,31]*/  jacp_[196] = kf_cAMPxRii__Rii_cAMP*(Rii_ConservedConst-(RiiP+RiiP_cAMP+Rii_cAMP+RiiP_CaN+RiiP_cAMP_CaN-C-AKAR4_C));
	/*[ 5,32]*/  jacp_[197] = kf_cAMPxRii__Rii_cAMP*(cAMP_ConservedConst-(RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN));
	/*[ 6,10]*/  jacp_[208] = (Rii_C_ConservedConst-(RiiP_C+RiiP_C_cAMP+C+Rii_C_cAMP+AKAR4_C))*(cAMP_ConservedConst-(RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN));
	/*[ 6,11]*/  jacp_[209] = -Rii_C_cAMP;
	/*[ 6,13]*/  jacp_[211] = Rii_cAMP*C;
	/*[ 6,14]*/  jacp_[212] = -Rii_C_cAMP;
	/*[ 6,15]*/  jacp_[213] = -Rii_C_cAMP;
	/*[ 6,30]*/  jacp_[228] = kf_Rii_CxcAMP__Rii_C_cAMP*(cAMP_ConservedConst-(RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN));
	/*[ 6,31]*/  jacp_[229] = kf_Rii_CxcAMP__Rii_C_cAMP*(Rii_C_ConservedConst-(RiiP_C+RiiP_C_cAMP+C+Rii_C_cAMP+AKAR4_C));
	/*[ 7,17]*/  jacp_[248] = ((b_AKAP*(1-b_AKAP))/kmON+gsl_pow_2(1-b_AKAP)/kmOFF)*RiiP*(CaN_ConservedConst-(RiiP_CaN+RiiP_cAMP_CaN))-(1-b_AKAP)*RiiP_CaN;
	/*[ 7,18]*/  jacp_[249] = ((b_AKAP*(1-b_AKAP))/kmON+gsl_pow_2(1-b_AKAP)/kmOFF)*RiiP*(CaN_ConservedConst-(RiiP_CaN+RiiP_cAMP_CaN))-(1-b_AKAP)*RiiP_CaN;
	/*[ 7,19]*/  jacp_[250] = (gsl_pow_2(b_AKAP)/kmON+((1-b_AKAP)*b_AKAP)/kmOFF)*RiiP*(CaN_ConservedConst-(RiiP_CaN+RiiP_cAMP_CaN))-b_AKAP*RiiP_CaN;
	/*[ 7,20]*/  jacp_[251] = (gsl_pow_2(b_AKAP)/kmON+((1-b_AKAP)*b_AKAP)/kmOFF)*RiiP*(CaN_ConservedConst-(RiiP_CaN+RiiP_cAMP_CaN))-b_AKAP*RiiP_CaN;
	/*[ 7,24]*/  jacp_[255] = (-(CaN_ConservedConst-(RiiP_CaN+RiiP_cAMP_CaN))*RiiP*(1-b_AKAP)*(b_AKAP*AKAPon_3+(1-b_AKAP)*AKAPoff_3+b_AKAP*AKAPon_1+(1-b_AKAP)*AKAPoff_1))/gsl_pow_2(kmOFF);
	/*[ 7,25]*/  jacp_[256] = (-(CaN_ConservedConst-(RiiP_CaN+RiiP_cAMP_CaN))*RiiP*b_AKAP*(b_AKAP*AKAPon_3+(1-b_AKAP)*AKAPoff_3+b_AKAP*AKAPon_1+(1-b_AKAP)*AKAPoff_1))/gsl_pow_2(kmON);
	/*[ 7,27]*/  jacp_[258] = ((b_AKAP*(AKAPon_3-AKAPoff_3+AKAPon_1-AKAPoff_1)+b_AKAP*AKAPon_3+(1-b_AKAP)*AKAPoff_3+b_AKAP*AKAPon_1+(1-b_AKAP)*AKAPoff_1)/kmON+((1-b_AKAP)*(AKAPon_3-AKAPoff_3+AKAPon_1-AKAPoff_1)-(b_AKAP*AKAPon_3+(1-b_AKAP)*AKAPoff_3+b_AKAP*AKAPon_1+(1-b_AKAP)*AKAPoff_1))/kmOFF)*RiiP*(CaN_ConservedConst-(RiiP_CaN+RiiP_cAMP_CaN))-(AKAPon_3-AKAPoff_3)*RiiP_CaN-(AKAPon_1-AKAPoff_1)*RiiP_CaN;
	/*[ 7,29]*/  jacp_[260] = ((b_AKAP*(b_AKAP*AKAPon_3+(1-b_AKAP)*AKAPoff_3+b_AKAP*AKAPon_1+(1-b_AKAP)*AKAPoff_1))/kmON+((1-b_AKAP)*(b_AKAP*AKAPon_3+(1-b_AKAP)*AKAPoff_3+b_AKAP*AKAPon_1+(1-b_AKAP)*AKAPoff_1))/kmOFF)*RiiP;
	/*[ 8,17]*/  jacp_[281] = ((b_AKAP*(1-b_AKAP))/kmON+gsl_pow_2(1-b_AKAP)/kmOFF)*(CaN_ConservedConst-(RiiP_CaN+RiiP_cAMP_CaN))*RiiP_cAMP-(1-b_AKAP)*RiiP_cAMP_CaN;
	/*[ 8,18]*/  jacp_[282] = ((b_AKAP*(1-b_AKAP))/kmON+gsl_pow_2(1-b_AKAP)/kmOFF)*(CaN_ConservedConst-(RiiP_CaN+RiiP_cAMP_CaN))*RiiP_cAMP-(1-b_AKAP)*RiiP_cAMP_CaN;
	/*[ 8,19]*/  jacp_[283] = (gsl_pow_2(b_AKAP)/kmON+((1-b_AKAP)*b_AKAP)/kmOFF)*(CaN_ConservedConst-(RiiP_CaN+RiiP_cAMP_CaN))*RiiP_cAMP-b_AKAP*RiiP_cAMP_CaN;
	/*[ 8,20]*/  jacp_[284] = (gsl_pow_2(b_AKAP)/kmON+((1-b_AKAP)*b_AKAP)/kmOFF)*(CaN_ConservedConst-(RiiP_CaN+RiiP_cAMP_CaN))*RiiP_cAMP-b_AKAP*RiiP_cAMP_CaN;
	/*[ 8,24]*/  jacp_[288] = (-RiiP_cAMP*(CaN_ConservedConst-(RiiP_CaN+RiiP_cAMP_CaN))*(1-b_AKAP)*(b_AKAP*AKAPon_3+(1-b_AKAP)*AKAPoff_3+b_AKAP*AKAPon_1+(1-b_AKAP)*AKAPoff_1))/gsl_pow_2(kmOFF);
	/*[ 8,25]*/  jacp_[289] = (-RiiP_cAMP*(CaN_ConservedConst-(RiiP_CaN+RiiP_cAMP_CaN))*b_AKAP*(b_AKAP*AKAPon_3+(1-b_AKAP)*AKAPoff_3+b_AKAP*AKAPon_1+(1-b_AKAP)*AKAPoff_1))/gsl_pow_2(kmON);
	/*[ 8,27]*/  jacp_[291] = ((b_AKAP*(AKAPon_3-AKAPoff_3+AKAPon_1-AKAPoff_1)+b_AKAP*AKAPon_3+(1-b_AKAP)*AKAPoff_3+b_AKAP*AKAPon_1+(1-b_AKAP)*AKAPoff_1)/kmON+((1-b_AKAP)*(AKAPon_3-AKAPoff_3+AKAPon_1-AKAPoff_1)-(b_AKAP*AKAPon_3+(1-b_AKAP)*AKAPoff_3+b_AKAP*AKAPon_1+(1-b_AKAP)*AKAPoff_1))/kmOFF)*(CaN_ConservedConst-(RiiP_CaN+RiiP_cAMP_CaN))*RiiP_cAMP-(AKAPon_3-AKAPoff_3)*RiiP_cAMP_CaN-(AKAPon_1-AKAPoff_1)*RiiP_cAMP_CaN;
	/*[ 8,29]*/  jacp_[293] = ((b_AKAP*(b_AKAP*AKAPon_3+(1-b_AKAP)*AKAPoff_3+b_AKAP*AKAPon_1+(1-b_AKAP)*AKAPoff_1))/kmON+((1-b_AKAP)*(b_AKAP*AKAPon_3+(1-b_AKAP)*AKAPoff_3+b_AKAP*AKAPon_1+(1-b_AKAP)*AKAPoff_1))/kmOFF)*RiiP_cAMP;
	/*[ 9,21]*/  jacp_[318] = C*(AKAR4_ConservedConst-(AKAR4_C+AKAR4p));
	/*[ 9,22]*/  jacp_[319] = -AKAR4_C;
	/*[ 9,23]*/  jacp_[320] = -AKAR4_C;
	/*[ 9,28]*/  jacp_[325] = kf_C_AKAR4*C;
	/*[10,23]*/  jacp_[353] = AKAR4_C;
	return GSL_SUCCESS;
}

/* Output Function (Observables)   */
int AKAP79_func(double t, const double y_[], double *func_, void *par){
	double *p_=par;
	if (!y_ || !func_) return 1;
/* 	constants   */
/* 	parameter values   */
	double kf_Rii_C__RiiP_C = p_[_kf_Rii_C__RiiP_C];                    /* [  0] */
	double kf_RiiP_CxcAMP__RiiP_C_cAMP = p_[_kf_RiiP_CxcAMP__RiiP_C_cAMP];       /* [  1] */
	double kf_RiiP_cAMPxC__RiiP_C_cAMP = p_[_kf_RiiP_cAMPxC__RiiP_C_cAMP];       /* [  2] */
	double kb_RiiP_cAMPxC__RiiP_C_cAMP = p_[_kb_RiiP_cAMPxC__RiiP_C_cAMP];       /* [  3] */
	double kb_RiiPXcAMP__RiiP_cAMP = p_[_kb_RiiPXcAMP__RiiP_cAMP];       /* [  4] */
	double kf_RiiPXcAMP__RiiP_cAMP = p_[_kf_RiiPXcAMP__RiiP_cAMP];       /* [  5] */
	double kf_RiiPxC__RiiP_C = p_[_kf_RiiPxC__RiiP_C];                  /* [  6] */
	double kb_RiiPxC__RiiP_C = p_[_kb_RiiPxC__RiiP_C];                  /* [  7] */
	double kf_cAMPxRii__Rii_cAMP = p_[_kf_cAMPxRii__Rii_cAMP];          /* [  8] */
	double kb_cAMPxRii__Rii_cAMP = p_[_kb_cAMPxRii__Rii_cAMP];          /* [  9] */
	double kf_Rii_CxcAMP__Rii_C_cAMP = p_[_kf_Rii_CxcAMP__Rii_C_cAMP];       /* [ 10] */
	double kb_Rii_CxcAMP__Rii_C_cAMP = p_[_kb_Rii_CxcAMP__Rii_C_cAMP];       /* [ 11] */
	double kf_RiixC__Rii_C = p_[_kf_RiixC__Rii_C];                      /* [ 12] */
	double kf_Rii_cAMPxC__Rii_C_cAMP = p_[_kf_Rii_cAMPxC__Rii_C_cAMP];       /* [ 13] */
	double kb_Rii_cAMPxC__Rii_C_cAMP = p_[_kb_Rii_cAMPxC__Rii_C_cAMP];       /* [ 14] */
	double kf_Rii_C_cAMP__RiiP_C_cAMP = p_[_kf_Rii_C_cAMP__RiiP_C_cAMP];       /* [ 15] */
	double kb_RiixC__Rii_C = p_[_kb_RiixC__Rii_C];                      /* [ 16] */
	double AKAPoff_1 = p_[_AKAPoff_1];                                  /* [ 17] */
	double AKAPoff_3 = p_[_AKAPoff_3];                                  /* [ 18] */
	double AKAPon_1 = p_[_AKAPon_1];                                    /* [ 19] */
	double AKAPon_3 = p_[_AKAPon_3];                                    /* [ 20] */
	double kf_C_AKAR4 = p_[_kf_C_AKAR4];                                /* [ 21] */
	double kb_C_AKAR4 = p_[_kb_C_AKAR4];                                /* [ 22] */
	double kcat_AKARp = p_[_kcat_AKARp];                                /* [ 23] */
	double kmOFF = p_[_kmOFF];                                          /* [ 24] */
	double kmON = p_[_kmON];                                            /* [ 25] */
	double KD_T = p_[_KD_T];                                            /* [ 26] */
	double b_AKAP = p_[_b_AKAP];                                        /* [ 27] */
	double AKAR4_ConservedConst = p_[_AKAR4_ConservedConst];            /* [ 28] */
	double CaN_ConservedConst = p_[_CaN_ConservedConst];                /* [ 29] */
	double Rii_C_ConservedConst = p_[_Rii_C_ConservedConst];            /* [ 30] */
	double cAMP_ConservedConst = p_[_cAMP_ConservedConst];              /* [ 31] */
	double Rii_ConservedConst = p_[_Rii_ConservedConst];                /* [ 32] */
/* 	state variables   */
	double RiiP = y_[_RiiP];                                            /* [  0] */
	double RiiP_cAMP = y_[_RiiP_cAMP];                                  /* [  1] */
	double RiiP_C = y_[_RiiP_C];                                        /* [  2] */
	double RiiP_C_cAMP = y_[_RiiP_C_cAMP];                              /* [  3] */
	double C = y_[_C];                                                  /* [  4] */
	double Rii_cAMP = y_[_Rii_cAMP];                                    /* [  5] */
	double Rii_C_cAMP = y_[_Rii_C_cAMP];                                /* [  6] */
	double RiiP_CaN = y_[_RiiP_CaN];                                    /* [  7] */
	double RiiP_cAMP_CaN = y_[_RiiP_cAMP_CaN];                          /* [  8] */
	double AKAR4_C = y_[_AKAR4_C];                                      /* [  9] */
	double AKAR4p = y_[_AKAR4p];                                        /* [ 10] */
/* 	expressions   */
	double kf_RiiP_cAMP_CaN__CaNXRii_cAMP = b_AKAP * AKAPon_1 + (1 - b_AKAP) * AKAPoff_1;
	double kb_RiiPxCaN__RiiP_CaN = b_AKAP*AKAPon_3  +  (1 - b_AKAP)* AKAPoff_3;
	double kf_RiiP_CaN__RiixCaN = b_AKAP * AKAPon_1 + (1 - b_AKAP) * AKAPoff_1;
	double kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN = b_AKAP*AKAPon_3  +  (1 - b_AKAP)* AKAPoff_3;
	double kf_RiiPxCaN__RiiP_CaN = b_AKAP * ((kb_RiiPxCaN__RiiP_CaN + kf_RiiP_cAMP_CaN__CaNXRii_cAMP)/kmON ) + (1 - b_AKAP) * (kb_RiiPxCaN__RiiP_CaN + kf_RiiP_cAMP_CaN__CaNXRii_cAMP) / kmOFF;
	double kf_CaNxRiiP_cAMP__RiiP_cAMP_CaN = b_AKAP * ((kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN + kf_RiiP_CaN__RiixCaN)/kmON) + (1 - b_AKAP) * (kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN + kf_RiiP_CaN__RiixCaN)/kmOFF;
	double kb_RiiP_CxcAMP__RiiP_C_cAMP = kf_RiiP_CxcAMP__RiiP_C_cAMP * KD_T;
	double AKAR4 = (AKAR4_ConservedConst - (AKAR4_C+AKAR4p));
	double CaN = (CaN_ConservedConst - (RiiP_CaN+RiiP_cAMP_CaN));
	double Rii_C = (Rii_C_ConservedConst - (RiiP_C+RiiP_C_cAMP+C+Rii_C_cAMP+AKAR4_C));
	double cAMP = (cAMP_ConservedConst - (RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN));
	double Rii = (Rii_ConservedConst - (RiiP+RiiP_cAMP+Rii_cAMP+RiiP_CaN+RiiP_cAMP_CaN-C-AKAR4_C));
	double reaction_51 = kf_Rii_C__RiiP_C*Rii_C;
	double reaction_14 = kf_RiiPxC__RiiP_C*RiiP*C - kb_RiiPxC__RiiP_C*RiiP_C;
	double reaction_12 = kf_RiiP_CxcAMP__RiiP_C_cAMP*RiiP_C*cAMP - kb_RiiP_CxcAMP__RiiP_C_cAMP*RiiP_C_cAMP;
	double reaction_43 = kf_RiiPXcAMP__RiiP_cAMP*cAMP*RiiP - kb_RiiPXcAMP__RiiP_cAMP*RiiP_cAMP;
	double reaction_23 = kf_RiiP_cAMPxC__RiiP_C_cAMP*RiiP_cAMP*C - kb_RiiP_cAMPxC__RiiP_C_cAMP*RiiP_C_cAMP;
	double reaction_78 = kf_cAMPxRii__Rii_cAMP*cAMP*Rii - kb_cAMPxRii__Rii_cAMP*Rii_cAMP;
	double reaction_56 = kf_Rii_CxcAMP__Rii_C_cAMP*Rii_C*cAMP - kb_Rii_CxcAMP__Rii_C_cAMP*Rii_C_cAMP;
	double reaction_76 = kf_Rii_cAMPxC__Rii_C_cAMP*Rii_cAMP*C - kb_Rii_cAMPxC__Rii_C_cAMP*Rii_C_cAMP;
	double reaction_62 = kf_Rii_C_cAMP__RiiP_C_cAMP*Rii_C_cAMP;
	double reaction_58 = kf_RiixC__Rii_C*Rii*C - kb_RiixC__Rii_C*Rii_C;
	double reaction_44 = kf_RiiPxCaN__RiiP_CaN*RiiP*CaN - kb_RiiPxCaN__RiiP_CaN*RiiP_CaN;
	double reaction_33 = kf_CaNxRiiP_cAMP__RiiP_cAMP_CaN*CaN*RiiP_cAMP - kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN*RiiP_cAMP_CaN;
	double reaction_48 = kf_RiiP_CaN__RiixCaN*RiiP_CaN;
	double reaction_37 = kf_RiiP_cAMP_CaN__CaNXRii_cAMP*RiiP_cAMP_CaN;
	double reaction_1 = kf_C_AKAR4*C*AKAR4 - kb_C_AKAR4*AKAR4_C;
	double reaction_2 = kcat_AKARp*AKAR4_C;
	func_[_AKAR4pOUT] = (AKAR4p*5)*71.67+100;
	return GSL_SUCCESS;
}

/* Output function Jacobian: dF(t,y;p)/dx   */
int AKAP79_funcJac(double t, const double y_[], double *funcJac_, void *par){
	double *p_=par;
	if (!y_ || !funcJac_) return 11;
/* 	constants   */
/* 	parameter values   */
	double kf_Rii_C__RiiP_C = p_[_kf_Rii_C__RiiP_C];                    /* [  0] */
	double kf_RiiP_CxcAMP__RiiP_C_cAMP = p_[_kf_RiiP_CxcAMP__RiiP_C_cAMP];       /* [  1] */
	double kf_RiiP_cAMPxC__RiiP_C_cAMP = p_[_kf_RiiP_cAMPxC__RiiP_C_cAMP];       /* [  2] */
	double kb_RiiP_cAMPxC__RiiP_C_cAMP = p_[_kb_RiiP_cAMPxC__RiiP_C_cAMP];       /* [  3] */
	double kb_RiiPXcAMP__RiiP_cAMP = p_[_kb_RiiPXcAMP__RiiP_cAMP];       /* [  4] */
	double kf_RiiPXcAMP__RiiP_cAMP = p_[_kf_RiiPXcAMP__RiiP_cAMP];       /* [  5] */
	double kf_RiiPxC__RiiP_C = p_[_kf_RiiPxC__RiiP_C];                  /* [  6] */
	double kb_RiiPxC__RiiP_C = p_[_kb_RiiPxC__RiiP_C];                  /* [  7] */
	double kf_cAMPxRii__Rii_cAMP = p_[_kf_cAMPxRii__Rii_cAMP];          /* [  8] */
	double kb_cAMPxRii__Rii_cAMP = p_[_kb_cAMPxRii__Rii_cAMP];          /* [  9] */
	double kf_Rii_CxcAMP__Rii_C_cAMP = p_[_kf_Rii_CxcAMP__Rii_C_cAMP];       /* [ 10] */
	double kb_Rii_CxcAMP__Rii_C_cAMP = p_[_kb_Rii_CxcAMP__Rii_C_cAMP];       /* [ 11] */
	double kf_RiixC__Rii_C = p_[_kf_RiixC__Rii_C];                      /* [ 12] */
	double kf_Rii_cAMPxC__Rii_C_cAMP = p_[_kf_Rii_cAMPxC__Rii_C_cAMP];       /* [ 13] */
	double kb_Rii_cAMPxC__Rii_C_cAMP = p_[_kb_Rii_cAMPxC__Rii_C_cAMP];       /* [ 14] */
	double kf_Rii_C_cAMP__RiiP_C_cAMP = p_[_kf_Rii_C_cAMP__RiiP_C_cAMP];       /* [ 15] */
	double kb_RiixC__Rii_C = p_[_kb_RiixC__Rii_C];                      /* [ 16] */
	double AKAPoff_1 = p_[_AKAPoff_1];                                  /* [ 17] */
	double AKAPoff_3 = p_[_AKAPoff_3];                                  /* [ 18] */
	double AKAPon_1 = p_[_AKAPon_1];                                    /* [ 19] */
	double AKAPon_3 = p_[_AKAPon_3];                                    /* [ 20] */
	double kf_C_AKAR4 = p_[_kf_C_AKAR4];                                /* [ 21] */
	double kb_C_AKAR4 = p_[_kb_C_AKAR4];                                /* [ 22] */
	double kcat_AKARp = p_[_kcat_AKARp];                                /* [ 23] */
	double kmOFF = p_[_kmOFF];                                          /* [ 24] */
	double kmON = p_[_kmON];                                            /* [ 25] */
	double KD_T = p_[_KD_T];                                            /* [ 26] */
	double b_AKAP = p_[_b_AKAP];                                        /* [ 27] */
	double AKAR4_ConservedConst = p_[_AKAR4_ConservedConst];            /* [ 28] */
	double CaN_ConservedConst = p_[_CaN_ConservedConst];                /* [ 29] */
	double Rii_C_ConservedConst = p_[_Rii_C_ConservedConst];            /* [ 30] */
	double cAMP_ConservedConst = p_[_cAMP_ConservedConst];              /* [ 31] */
	double Rii_ConservedConst = p_[_Rii_ConservedConst];                /* [ 32] */
/* 	state variables   */
	double RiiP = y_[_RiiP];                                            /* [  0] */
	double RiiP_cAMP = y_[_RiiP_cAMP];                                  /* [  1] */
	double RiiP_C = y_[_RiiP_C];                                        /* [  2] */
	double RiiP_C_cAMP = y_[_RiiP_C_cAMP];                              /* [  3] */
	double C = y_[_C];                                                  /* [  4] */
	double Rii_cAMP = y_[_Rii_cAMP];                                    /* [  5] */
	double Rii_C_cAMP = y_[_Rii_C_cAMP];                                /* [  6] */
	double RiiP_CaN = y_[_RiiP_CaN];                                    /* [  7] */
	double RiiP_cAMP_CaN = y_[_RiiP_cAMP_CaN];                          /* [  8] */
	double AKAR4_C = y_[_AKAR4_C];                                      /* [  9] */
	double AKAR4p = y_[_AKAR4p];                                        /* [ 10] */
/* 	expressions   */
	double kf_RiiP_cAMP_CaN__CaNXRii_cAMP = b_AKAP * AKAPon_1 + (1 - b_AKAP) * AKAPoff_1;
	double kb_RiiPxCaN__RiiP_CaN = b_AKAP*AKAPon_3  +  (1 - b_AKAP)* AKAPoff_3;
	double kf_RiiP_CaN__RiixCaN = b_AKAP * AKAPon_1 + (1 - b_AKAP) * AKAPoff_1;
	double kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN = b_AKAP*AKAPon_3  +  (1 - b_AKAP)* AKAPoff_3;
	double kf_RiiPxCaN__RiiP_CaN = b_AKAP * ((kb_RiiPxCaN__RiiP_CaN + kf_RiiP_cAMP_CaN__CaNXRii_cAMP)/kmON ) + (1 - b_AKAP) * (kb_RiiPxCaN__RiiP_CaN + kf_RiiP_cAMP_CaN__CaNXRii_cAMP) / kmOFF;
	double kf_CaNxRiiP_cAMP__RiiP_cAMP_CaN = b_AKAP * ((kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN + kf_RiiP_CaN__RiixCaN)/kmON) + (1 - b_AKAP) * (kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN + kf_RiiP_CaN__RiixCaN)/kmOFF;
	double kb_RiiP_CxcAMP__RiiP_C_cAMP = kf_RiiP_CxcAMP__RiiP_C_cAMP * KD_T;
	double AKAR4 = (AKAR4_ConservedConst - (AKAR4_C+AKAR4p));
	double CaN = (CaN_ConservedConst - (RiiP_CaN+RiiP_cAMP_CaN));
	double Rii_C = (Rii_C_ConservedConst - (RiiP_C+RiiP_C_cAMP+C+Rii_C_cAMP+AKAR4_C));
	double cAMP = (cAMP_ConservedConst - (RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN));
	double Rii = (Rii_ConservedConst - (RiiP+RiiP_cAMP+Rii_cAMP+RiiP_CaN+RiiP_cAMP_CaN-C-AKAR4_C));
	double reaction_51 = kf_Rii_C__RiiP_C*Rii_C;
	double reaction_14 = kf_RiiPxC__RiiP_C*RiiP*C - kb_RiiPxC__RiiP_C*RiiP_C;
	double reaction_12 = kf_RiiP_CxcAMP__RiiP_C_cAMP*RiiP_C*cAMP - kb_RiiP_CxcAMP__RiiP_C_cAMP*RiiP_C_cAMP;
	double reaction_43 = kf_RiiPXcAMP__RiiP_cAMP*cAMP*RiiP - kb_RiiPXcAMP__RiiP_cAMP*RiiP_cAMP;
	double reaction_23 = kf_RiiP_cAMPxC__RiiP_C_cAMP*RiiP_cAMP*C - kb_RiiP_cAMPxC__RiiP_C_cAMP*RiiP_C_cAMP;
	double reaction_78 = kf_cAMPxRii__Rii_cAMP*cAMP*Rii - kb_cAMPxRii__Rii_cAMP*Rii_cAMP;
	double reaction_56 = kf_Rii_CxcAMP__Rii_C_cAMP*Rii_C*cAMP - kb_Rii_CxcAMP__Rii_C_cAMP*Rii_C_cAMP;
	double reaction_76 = kf_Rii_cAMPxC__Rii_C_cAMP*Rii_cAMP*C - kb_Rii_cAMPxC__Rii_C_cAMP*Rii_C_cAMP;
	double reaction_62 = kf_Rii_C_cAMP__RiiP_C_cAMP*Rii_C_cAMP;
	double reaction_58 = kf_RiixC__Rii_C*Rii*C - kb_RiixC__Rii_C*Rii_C;
	double reaction_44 = kf_RiiPxCaN__RiiP_CaN*RiiP*CaN - kb_RiiPxCaN__RiiP_CaN*RiiP_CaN;
	double reaction_33 = kf_CaNxRiiP_cAMP__RiiP_cAMP_CaN*CaN*RiiP_cAMP - kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN*RiiP_cAMP_CaN;
	double reaction_48 = kf_RiiP_CaN__RiixCaN*RiiP_CaN;
	double reaction_37 = kf_RiiP_cAMP_CaN__CaNXRii_cAMP*RiiP_cAMP_CaN;
	double reaction_1 = kf_C_AKAR4*C*AKAR4 - kb_C_AKAR4*AKAR4_C;
	double reaction_2 = kcat_AKARp*AKAR4_C;
	memset(funcJac_,0,sizeof(double)*11); /* initialize with 0.0 */
	/*[ 0,10]*/  funcJac_[10] = 358.35;
	return GSL_SUCCESS;
}

/* Output function parameter Jacobian: dF(t,y;p)/dp   */
int AKAP79_funcJacp(double t, const double y_[], double *funcJacp_, void *par){
	double *p_=par;
	if (!y_ || !funcJacp_) return 33;
/* 	constants   */
/* 	parameter values   */
	double kf_Rii_C__RiiP_C = p_[_kf_Rii_C__RiiP_C];                    /* [  0] */
	double kf_RiiP_CxcAMP__RiiP_C_cAMP = p_[_kf_RiiP_CxcAMP__RiiP_C_cAMP];       /* [  1] */
	double kf_RiiP_cAMPxC__RiiP_C_cAMP = p_[_kf_RiiP_cAMPxC__RiiP_C_cAMP];       /* [  2] */
	double kb_RiiP_cAMPxC__RiiP_C_cAMP = p_[_kb_RiiP_cAMPxC__RiiP_C_cAMP];       /* [  3] */
	double kb_RiiPXcAMP__RiiP_cAMP = p_[_kb_RiiPXcAMP__RiiP_cAMP];       /* [  4] */
	double kf_RiiPXcAMP__RiiP_cAMP = p_[_kf_RiiPXcAMP__RiiP_cAMP];       /* [  5] */
	double kf_RiiPxC__RiiP_C = p_[_kf_RiiPxC__RiiP_C];                  /* [  6] */
	double kb_RiiPxC__RiiP_C = p_[_kb_RiiPxC__RiiP_C];                  /* [  7] */
	double kf_cAMPxRii__Rii_cAMP = p_[_kf_cAMPxRii__Rii_cAMP];          /* [  8] */
	double kb_cAMPxRii__Rii_cAMP = p_[_kb_cAMPxRii__Rii_cAMP];          /* [  9] */
	double kf_Rii_CxcAMP__Rii_C_cAMP = p_[_kf_Rii_CxcAMP__Rii_C_cAMP];       /* [ 10] */
	double kb_Rii_CxcAMP__Rii_C_cAMP = p_[_kb_Rii_CxcAMP__Rii_C_cAMP];       /* [ 11] */
	double kf_RiixC__Rii_C = p_[_kf_RiixC__Rii_C];                      /* [ 12] */
	double kf_Rii_cAMPxC__Rii_C_cAMP = p_[_kf_Rii_cAMPxC__Rii_C_cAMP];       /* [ 13] */
	double kb_Rii_cAMPxC__Rii_C_cAMP = p_[_kb_Rii_cAMPxC__Rii_C_cAMP];       /* [ 14] */
	double kf_Rii_C_cAMP__RiiP_C_cAMP = p_[_kf_Rii_C_cAMP__RiiP_C_cAMP];       /* [ 15] */
	double kb_RiixC__Rii_C = p_[_kb_RiixC__Rii_C];                      /* [ 16] */
	double AKAPoff_1 = p_[_AKAPoff_1];                                  /* [ 17] */
	double AKAPoff_3 = p_[_AKAPoff_3];                                  /* [ 18] */
	double AKAPon_1 = p_[_AKAPon_1];                                    /* [ 19] */
	double AKAPon_3 = p_[_AKAPon_3];                                    /* [ 20] */
	double kf_C_AKAR4 = p_[_kf_C_AKAR4];                                /* [ 21] */
	double kb_C_AKAR4 = p_[_kb_C_AKAR4];                                /* [ 22] */
	double kcat_AKARp = p_[_kcat_AKARp];                                /* [ 23] */
	double kmOFF = p_[_kmOFF];                                          /* [ 24] */
	double kmON = p_[_kmON];                                            /* [ 25] */
	double KD_T = p_[_KD_T];                                            /* [ 26] */
	double b_AKAP = p_[_b_AKAP];                                        /* [ 27] */
	double AKAR4_ConservedConst = p_[_AKAR4_ConservedConst];            /* [ 28] */
	double CaN_ConservedConst = p_[_CaN_ConservedConst];                /* [ 29] */
	double Rii_C_ConservedConst = p_[_Rii_C_ConservedConst];            /* [ 30] */
	double cAMP_ConservedConst = p_[_cAMP_ConservedConst];              /* [ 31] */
	double Rii_ConservedConst = p_[_Rii_ConservedConst];                /* [ 32] */
/* 	state variables   */
	double RiiP = y_[_RiiP];                                            /* [  0] */
	double RiiP_cAMP = y_[_RiiP_cAMP];                                  /* [  1] */
	double RiiP_C = y_[_RiiP_C];                                        /* [  2] */
	double RiiP_C_cAMP = y_[_RiiP_C_cAMP];                              /* [  3] */
	double C = y_[_C];                                                  /* [  4] */
	double Rii_cAMP = y_[_Rii_cAMP];                                    /* [  5] */
	double Rii_C_cAMP = y_[_Rii_C_cAMP];                                /* [  6] */
	double RiiP_CaN = y_[_RiiP_CaN];                                    /* [  7] */
	double RiiP_cAMP_CaN = y_[_RiiP_cAMP_CaN];                          /* [  8] */
	double AKAR4_C = y_[_AKAR4_C];                                      /* [  9] */
	double AKAR4p = y_[_AKAR4p];                                        /* [ 10] */
/* 	expressions   */
	double kf_RiiP_cAMP_CaN__CaNXRii_cAMP = b_AKAP * AKAPon_1 + (1 - b_AKAP) * AKAPoff_1;
	double kb_RiiPxCaN__RiiP_CaN = b_AKAP*AKAPon_3  +  (1 - b_AKAP)* AKAPoff_3;
	double kf_RiiP_CaN__RiixCaN = b_AKAP * AKAPon_1 + (1 - b_AKAP) * AKAPoff_1;
	double kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN = b_AKAP*AKAPon_3  +  (1 - b_AKAP)* AKAPoff_3;
	double kf_RiiPxCaN__RiiP_CaN = b_AKAP * ((kb_RiiPxCaN__RiiP_CaN + kf_RiiP_cAMP_CaN__CaNXRii_cAMP)/kmON ) + (1 - b_AKAP) * (kb_RiiPxCaN__RiiP_CaN + kf_RiiP_cAMP_CaN__CaNXRii_cAMP) / kmOFF;
	double kf_CaNxRiiP_cAMP__RiiP_cAMP_CaN = b_AKAP * ((kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN + kf_RiiP_CaN__RiixCaN)/kmON) + (1 - b_AKAP) * (kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN + kf_RiiP_CaN__RiixCaN)/kmOFF;
	double kb_RiiP_CxcAMP__RiiP_C_cAMP = kf_RiiP_CxcAMP__RiiP_C_cAMP * KD_T;
	double AKAR4 = (AKAR4_ConservedConst - (AKAR4_C+AKAR4p));
	double CaN = (CaN_ConservedConst - (RiiP_CaN+RiiP_cAMP_CaN));
	double Rii_C = (Rii_C_ConservedConst - (RiiP_C+RiiP_C_cAMP+C+Rii_C_cAMP+AKAR4_C));
	double cAMP = (cAMP_ConservedConst - (RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN));
	double Rii = (Rii_ConservedConst - (RiiP+RiiP_cAMP+Rii_cAMP+RiiP_CaN+RiiP_cAMP_CaN-C-AKAR4_C));
	double reaction_51 = kf_Rii_C__RiiP_C*Rii_C;
	double reaction_14 = kf_RiiPxC__RiiP_C*RiiP*C - kb_RiiPxC__RiiP_C*RiiP_C;
	double reaction_12 = kf_RiiP_CxcAMP__RiiP_C_cAMP*RiiP_C*cAMP - kb_RiiP_CxcAMP__RiiP_C_cAMP*RiiP_C_cAMP;
	double reaction_43 = kf_RiiPXcAMP__RiiP_cAMP*cAMP*RiiP - kb_RiiPXcAMP__RiiP_cAMP*RiiP_cAMP;
	double reaction_23 = kf_RiiP_cAMPxC__RiiP_C_cAMP*RiiP_cAMP*C - kb_RiiP_cAMPxC__RiiP_C_cAMP*RiiP_C_cAMP;
	double reaction_78 = kf_cAMPxRii__Rii_cAMP*cAMP*Rii - kb_cAMPxRii__Rii_cAMP*Rii_cAMP;
	double reaction_56 = kf_Rii_CxcAMP__Rii_C_cAMP*Rii_C*cAMP - kb_Rii_CxcAMP__Rii_C_cAMP*Rii_C_cAMP;
	double reaction_76 = kf_Rii_cAMPxC__Rii_C_cAMP*Rii_cAMP*C - kb_Rii_cAMPxC__Rii_C_cAMP*Rii_C_cAMP;
	double reaction_62 = kf_Rii_C_cAMP__RiiP_C_cAMP*Rii_C_cAMP;
	double reaction_58 = kf_RiixC__Rii_C*Rii*C - kb_RiixC__Rii_C*Rii_C;
	double reaction_44 = kf_RiiPxCaN__RiiP_CaN*RiiP*CaN - kb_RiiPxCaN__RiiP_CaN*RiiP_CaN;
	double reaction_33 = kf_CaNxRiiP_cAMP__RiiP_cAMP_CaN*CaN*RiiP_cAMP - kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN*RiiP_cAMP_CaN;
	double reaction_48 = kf_RiiP_CaN__RiixCaN*RiiP_CaN;
	double reaction_37 = kf_RiiP_cAMP_CaN__CaNXRii_cAMP*RiiP_cAMP_CaN;
	double reaction_1 = kf_C_AKAR4*C*AKAR4 - kb_C_AKAR4*AKAR4_C;
	double reaction_2 = kcat_AKARp*AKAR4_C;
	memset(funcJacp_,0,sizeof(double)*33); /* initialize with 0.0 */
	return GSL_SUCCESS;
}

int AKAP79_default(double t, double *p_){
	if (!p_) return numParam;
/* 	constants   */
	memset(p_,0,sizeof(double)*33); /* initialize with 0.0 */
	p_[_kf_Rii_C__RiiP_C] = 33;
	p_[_kf_RiiP_CxcAMP__RiiP_C_cAMP] = 0.496;
	p_[_kf_RiiP_cAMPxC__RiiP_C_cAMP] = 0.00545;
	p_[_kb_RiiP_cAMPxC__RiiP_C_cAMP] = 0.0156;
	p_[_kb_RiiPXcAMP__RiiP_cAMP] = 0.0016;
	p_[_kf_RiiPXcAMP__RiiP_cAMP] = 0.015;
	p_[_kf_RiiPxC__RiiP_C] = 0.038;
	p_[_kb_RiiPxC__RiiP_C] = 0.0026;
	p_[_kf_cAMPxRii__Rii_cAMP] = 0.015;
	p_[_kb_cAMPxRii__Rii_cAMP] = 0.0016;
	p_[_kf_Rii_CxcAMP__Rii_C_cAMP] = 0.496;
	p_[_kb_Rii_CxcAMP__Rii_C_cAMP] = 1.413;
	p_[_kf_RiixC__Rii_C] = 2.1;
	p_[_kf_Rii_cAMPxC__Rii_C_cAMP] = 0.2984;
	p_[_kb_Rii_cAMPxC__Rii_C_cAMP] = 0.018;
	p_[_kf_Rii_C_cAMP__RiiP_C_cAMP] = 33;
	p_[_kb_RiixC__Rii_C] = 3e-04;
	p_[_AKAPoff_1] = 2.6;
	p_[_AKAPoff_3] = 20;
	p_[_AKAPon_1] = 0.45;
	p_[_AKAPon_3] = 2;
	p_[_kf_C_AKAR4] = 0.018;
	p_[_kb_C_AKAR4] = 0.106;
	p_[_kcat_AKARp] = 10.2;
	p_[_kmOFF] = 100;
	p_[_kmON] = 1;
	p_[_KD_T] = 0.7;
	p_[_AKAR4_ConservedConst] = 0.2;
	p_[_CaN_ConservedConst] = 1.5;
	p_[_Rii_C_ConservedConst] = 0.63;
	p_[_Rii_ConservedConst] = 6.3;
	return GSL_SUCCESS;
}

int AKAP79_init(double t, double *y_, void *par){
	double *p_=par;
	if (!y_ || !y_) return 11;
/* 	constants   */
/* 	parameter values   */
	double kf_Rii_C__RiiP_C = p_[_kf_Rii_C__RiiP_C];                    /* [  0] */
	double kf_RiiP_CxcAMP__RiiP_C_cAMP = p_[_kf_RiiP_CxcAMP__RiiP_C_cAMP];       /* [  1] */
	double kf_RiiP_cAMPxC__RiiP_C_cAMP = p_[_kf_RiiP_cAMPxC__RiiP_C_cAMP];       /* [  2] */
	double kb_RiiP_cAMPxC__RiiP_C_cAMP = p_[_kb_RiiP_cAMPxC__RiiP_C_cAMP];       /* [  3] */
	double kb_RiiPXcAMP__RiiP_cAMP = p_[_kb_RiiPXcAMP__RiiP_cAMP];       /* [  4] */
	double kf_RiiPXcAMP__RiiP_cAMP = p_[_kf_RiiPXcAMP__RiiP_cAMP];       /* [  5] */
	double kf_RiiPxC__RiiP_C = p_[_kf_RiiPxC__RiiP_C];                  /* [  6] */
	double kb_RiiPxC__RiiP_C = p_[_kb_RiiPxC__RiiP_C];                  /* [  7] */
	double kf_cAMPxRii__Rii_cAMP = p_[_kf_cAMPxRii__Rii_cAMP];          /* [  8] */
	double kb_cAMPxRii__Rii_cAMP = p_[_kb_cAMPxRii__Rii_cAMP];          /* [  9] */
	double kf_Rii_CxcAMP__Rii_C_cAMP = p_[_kf_Rii_CxcAMP__Rii_C_cAMP];       /* [ 10] */
	double kb_Rii_CxcAMP__Rii_C_cAMP = p_[_kb_Rii_CxcAMP__Rii_C_cAMP];       /* [ 11] */
	double kf_RiixC__Rii_C = p_[_kf_RiixC__Rii_C];                      /* [ 12] */
	double kf_Rii_cAMPxC__Rii_C_cAMP = p_[_kf_Rii_cAMPxC__Rii_C_cAMP];       /* [ 13] */
	double kb_Rii_cAMPxC__Rii_C_cAMP = p_[_kb_Rii_cAMPxC__Rii_C_cAMP];       /* [ 14] */
	double kf_Rii_C_cAMP__RiiP_C_cAMP = p_[_kf_Rii_C_cAMP__RiiP_C_cAMP];       /* [ 15] */
	double kb_RiixC__Rii_C = p_[_kb_RiixC__Rii_C];                      /* [ 16] */
	double AKAPoff_1 = p_[_AKAPoff_1];                                  /* [ 17] */
	double AKAPoff_3 = p_[_AKAPoff_3];                                  /* [ 18] */
	double AKAPon_1 = p_[_AKAPon_1];                                    /* [ 19] */
	double AKAPon_3 = p_[_AKAPon_3];                                    /* [ 20] */
	double kf_C_AKAR4 = p_[_kf_C_AKAR4];                                /* [ 21] */
	double kb_C_AKAR4 = p_[_kb_C_AKAR4];                                /* [ 22] */
	double kcat_AKARp = p_[_kcat_AKARp];                                /* [ 23] */
	double kmOFF = p_[_kmOFF];                                          /* [ 24] */
	double kmON = p_[_kmON];                                            /* [ 25] */
	double KD_T = p_[_KD_T];                                            /* [ 26] */
	double b_AKAP = p_[_b_AKAP];                                        /* [ 27] */
	double AKAR4_ConservedConst = p_[_AKAR4_ConservedConst];            /* [ 28] */
	double CaN_ConservedConst = p_[_CaN_ConservedConst];                /* [ 29] */
	double Rii_C_ConservedConst = p_[_Rii_C_ConservedConst];            /* [ 30] */
	double cAMP_ConservedConst = p_[_cAMP_ConservedConst];              /* [ 31] */
	double Rii_ConservedConst = p_[_Rii_ConservedConst];                /* [ 32] */
	memset(y_,0,sizeof(double)*11); /* initialize with 0.0 */
	return GSL_SUCCESS;
}

