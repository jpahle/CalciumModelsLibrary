// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// sim_ano
DataFrame sim_ano(DataFrame user_input_df, List user_sim_params, List user_model_params);
RcppExport SEXP _CalciumModelsLibrary_sim_ano(SEXP user_input_dfSEXP, SEXP user_sim_paramsSEXP, SEXP user_model_paramsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame >::type user_input_df(user_input_dfSEXP);
    Rcpp::traits::input_parameter< List >::type user_sim_params(user_sim_paramsSEXP);
    Rcpp::traits::input_parameter< List >::type user_model_params(user_model_paramsSEXP);
    rcpp_result_gen = Rcpp::wrap(sim_ano(user_input_df, user_sim_params, user_model_params));
    return rcpp_result_gen;
END_RCPP
}
// sim_calcineurin
DataFrame sim_calcineurin(DataFrame user_input_df, List user_sim_params, List user_model_params);
RcppExport SEXP _CalciumModelsLibrary_sim_calcineurin(SEXP user_input_dfSEXP, SEXP user_sim_paramsSEXP, SEXP user_model_paramsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame >::type user_input_df(user_input_dfSEXP);
    Rcpp::traits::input_parameter< List >::type user_sim_params(user_sim_paramsSEXP);
    Rcpp::traits::input_parameter< List >::type user_model_params(user_model_paramsSEXP);
    rcpp_result_gen = Rcpp::wrap(sim_calcineurin(user_input_df, user_sim_params, user_model_params));
    return rcpp_result_gen;
END_RCPP
}
// sim_calmodulin
DataFrame sim_calmodulin(DataFrame user_input_df, List user_sim_params, List user_model_params);
RcppExport SEXP _CalciumModelsLibrary_sim_calmodulin(SEXP user_input_dfSEXP, SEXP user_sim_paramsSEXP, SEXP user_model_paramsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame >::type user_input_df(user_input_dfSEXP);
    Rcpp::traits::input_parameter< List >::type user_sim_params(user_sim_paramsSEXP);
    Rcpp::traits::input_parameter< List >::type user_model_params(user_model_paramsSEXP);
    rcpp_result_gen = Rcpp::wrap(sim_calmodulin(user_input_df, user_sim_params, user_model_params));
    return rcpp_result_gen;
END_RCPP
}
// sim_camkii
DataFrame sim_camkii(DataFrame user_input_df, List user_sim_params, List user_model_params);
RcppExport SEXP _CalciumModelsLibrary_sim_camkii(SEXP user_input_dfSEXP, SEXP user_sim_paramsSEXP, SEXP user_model_paramsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame >::type user_input_df(user_input_dfSEXP);
    Rcpp::traits::input_parameter< List >::type user_sim_params(user_sim_paramsSEXP);
    Rcpp::traits::input_parameter< List >::type user_model_params(user_model_paramsSEXP);
    rcpp_result_gen = Rcpp::wrap(sim_camkii(user_input_df, user_sim_params, user_model_params));
    return rcpp_result_gen;
END_RCPP
}
// sim_glycphos
DataFrame sim_glycphos(DataFrame user_input_df, List user_sim_params, List user_model_params);
RcppExport SEXP _CalciumModelsLibrary_sim_glycphos(SEXP user_input_dfSEXP, SEXP user_sim_paramsSEXP, SEXP user_model_paramsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame >::type user_input_df(user_input_dfSEXP);
    Rcpp::traits::input_parameter< List >::type user_sim_params(user_sim_paramsSEXP);
    Rcpp::traits::input_parameter< List >::type user_model_params(user_model_paramsSEXP);
    rcpp_result_gen = Rcpp::wrap(sim_glycphos(user_input_df, user_sim_params, user_model_params));
    return rcpp_result_gen;
END_RCPP
}
// sim_pkc
DataFrame sim_pkc(DataFrame user_input_df, List user_sim_params, List user_model_params);
RcppExport SEXP _CalciumModelsLibrary_sim_pkc(SEXP user_input_dfSEXP, SEXP user_sim_paramsSEXP, SEXP user_model_paramsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame >::type user_input_df(user_input_dfSEXP);
    Rcpp::traits::input_parameter< List >::type user_sim_params(user_sim_paramsSEXP);
    Rcpp::traits::input_parameter< List >::type user_model_params(user_model_paramsSEXP);
    rcpp_result_gen = Rcpp::wrap(sim_pkc(user_input_df, user_sim_params, user_model_params));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_CalciumModelsLibrary_sim_ano", (DL_FUNC) &_CalciumModelsLibrary_sim_ano, 3},
    {"_CalciumModelsLibrary_sim_calcineurin", (DL_FUNC) &_CalciumModelsLibrary_sim_calcineurin, 3},
    {"_CalciumModelsLibrary_sim_calmodulin", (DL_FUNC) &_CalciumModelsLibrary_sim_calmodulin, 3},
    {"_CalciumModelsLibrary_sim_camkii", (DL_FUNC) &_CalciumModelsLibrary_sim_camkii, 3},
    {"_CalciumModelsLibrary_sim_glycphos", (DL_FUNC) &_CalciumModelsLibrary_sim_glycphos, 3},
    {"_CalciumModelsLibrary_sim_pkc", (DL_FUNC) &_CalciumModelsLibrary_sim_pkc, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_CalciumModelsLibrary(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
