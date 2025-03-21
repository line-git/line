#include <iostream>
#include <cstring>
#include "mpc.h"

using namespace std;

#include "global_vars.h"


void mpfr_set_str_rat(
  mpfr_t *mpfr_num,
  const char *rat_str
) {
  const char *slash_pos = nullptr;
  for (const char *p = rat_str; *p != '\0'; ++p) {
    if (*p == '/') {
      slash_pos = p;
      break;
    }
  }

  if (slash_pos) {
    size_t num_len = slash_pos - rat_str;
    char *num_str = (char *)malloc(num_len + 1);
    strncpy(num_str, rat_str, num_len);
    num_str[num_len] = '\0';

    mpfr_set_str(*mpfr_num, num_str, 10, MPFR_RNDN);
    free(num_str);

    mpfr_t den;
    mpfr_init2(den, wp2);
    mpfr_set_str(den, slash_pos + 1, 10, MPFR_RNDN);
    mpfr_div(*mpfr_num, *mpfr_num, den, MPFR_RNDN);
    mpfr_clear(den);
  } else {
    mpfr_set_str(*mpfr_num, rat_str, 10, MPFR_RNDN);
  }
}


void mpc_set_str_rat(
  mpc_t *mpc_num,
  char *rat_str
) {
  char *inner_str;
  int i = 0, count = 0, cmplx = 0, max_len = 10000;
  // cout << rat_str << endl;

  if (rat_str[i] == '(') {
    // if it starts with '(', it's a complex number:
    // skip '(' and read until white space is found
    cmplx = 1;
    i++;
    inner_str = (char*) malloc((strlen(rat_str)+1)*sizeof(char));
    // cout << "looking for ' '" << endl;
    while (rat_str[i] != ' ') {
      // cout << "rat_str[" << i << "] = " << rat_str[i] << endl;
      inner_str[count++] = rat_str[i++];
      if (i == max_len) {
        perror("error while parsing complex number: white space not found or real part too long");
        exit(1);
      }
    }
    inner_str[count] = '\0';
    // cout << "inner str = " << inner_str << endl;
    i++;
  } else {
    // otherwise, no need to extract real and imaginary parts
    cmplx = 0;
    inner_str = rat_str;
  }

  if (cmplx) {
    // convert real part
    mpfr_set_str_rat(&mpc_realref(*mpc_num), inner_str);
    // read imaginary part
    count = 0;
    // cout << "looking for ')'" << endl;
    while (rat_str[i] != ')') {
      // cout << "rat_str[" << i << "] = " << rat_str[i] << endl;
      inner_str[count++] = rat_str[i++];
      if (count == max_len) {
        perror("error while parsing complex number: right delimiter not found or imaginary part too long");
        exit(1);
      }
    }
    inner_str[count] = '\0';
    // cout << "inner str = " << inner_str << endl;
    // convert imaginary part
    mpfr_set_str_rat(&mpc_imagref(*mpc_num), inner_str);  
  } else {
    // convert real part
    mpfr_set_str_rat(&mpc_realref(*mpc_num), inner_str);
    // set imaginary part to zero
    mpfr_set_ui(mpc_imagref(*mpc_num), 0, MPFR_RNDN);
  }

}

