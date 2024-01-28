#pragma once
#include "type.h"

namespace rtfmm
{
/**
 * @brief Compute Pk=Qk*Gk for 8x8 child-interactions in 1 parent-interaction pair(IN0 -> OUT0),
 * every time compute a 2x2 child interaction
*/
inline void c2c_8x8x1(real*& M_, real*& IN0, real*& OUT0)
{
    real out_reg000, out_reg001, out_reg010, out_reg011;
    real  in_reg000,  in_reg001,  in_reg010,  in_reg011;
    real   m_reg000,   m_reg001,   m_reg010,   m_reg011;
    real   m_reg100,   m_reg101,   m_reg110,   m_reg111;
    for(int i1=0;i1<8;i1+=2){
      real* IN0_=IN0;
      out_reg000=OUT0[ 0]; out_reg001=OUT0[ 1];
      out_reg010=OUT0[ 2]; out_reg011=OUT0[ 3];
      for(int i2=0;i2<8;i2+=2){
        m_reg000=M_[ 0]; m_reg001=M_[ 1];
        m_reg010=M_[ 2]; m_reg011=M_[ 3];
        m_reg100=M_[16]; m_reg101=M_[17];
        m_reg110=M_[18]; m_reg111=M_[19];

        in_reg000=IN0_[0]; in_reg001=IN0_[1];
        in_reg010=IN0_[2]; in_reg011=IN0_[3];

        out_reg000 += m_reg000*in_reg000 - m_reg001*in_reg001;
        out_reg001 += m_reg000*in_reg001 + m_reg001*in_reg000;
        out_reg010 += m_reg010*in_reg000 - m_reg011*in_reg001;
        out_reg011 += m_reg010*in_reg001 + m_reg011*in_reg000;

        out_reg000 += m_reg100*in_reg010 - m_reg101*in_reg011;
        out_reg001 += m_reg100*in_reg011 + m_reg101*in_reg010;
        out_reg010 += m_reg110*in_reg010 - m_reg111*in_reg011;
        out_reg011 += m_reg110*in_reg011 + m_reg111*in_reg010;

        M_+=32; // Jump to (column+2).
        IN0_+=4;
      }
      OUT0[ 0]=out_reg000; OUT0[ 1]=out_reg001;
      OUT0[ 2]=out_reg010; OUT0[ 3]=out_reg011;
      M_+=4-64*2; // Jump back to first column (row+2).
      OUT0+=4;
    }
}

/**
 * @brief Compute Pk=Qk*Gk for 8x8 child-interactions in 1 parent-interaction pair(IN0 -> OUT0),
 * every time compute a 1x1 child interaction
*/
inline void c2c_8x8x1_naive(real*& M_, real*& IN, real*& OUT)
{
    real pkr, pki;
    real qkr, qki;
    real gkr, gki;
    for(int idx_out = 0; idx_out < 8; idx_out++)
    {
        real* IN0_ = IN;
        pkr = OUT[0]; 
        pki = OUT[1];
        for(int idx_in = 0; idx_in < 8; idx_in++)
        {
            gkr = M_[ 0]; 
            gki = M_[ 1];
            qkr = IN0_[0]; 
            qki = IN0_[1];
            pkr += gkr * qkr - gki * qki;
            pki += gkr * qki + gki * qkr; 
            M_ += 2 * 8; // move to next row
            IN0_ += 2; // move to next in complex
        }
        OUT[0] = pkr; 
        OUT[1] = pki;
        M_ += -2 * 8 * 8 + 2; // move to next colum
        OUT += 2; // move to next out complex
    }
}

/**
 * @brief Compute Pk=Qk*Gk for 1x27 child-interactions in 1 parent-interaction pair(IN0 -> OUT0),
 * every time compute a 1x1 child interaction
*/
inline void c2c_1x27x1_naive(real*& M_, real*& IN, real*& OUT)
{
    real qkr = IN[0];
    real qki = IN[1];
    real pkr = OUT[0]; 
    real pki = OUT[1];
    real gkr, gki;
    for(int i = 0; i < 27; i++)
    {
        gkr = M_[0]; 
        gki = M_[1];
        pkr += gkr * qkr - gki * qki;
        pki += gkr * qki + gki * qkr; 
        M_ += 2;
    }
    OUT[0] = pkr; 
    OUT[1] = pki;
}

/**
 * @brief Compute Pk=Qk*Gk for 1x27 child-interactions in 1 parent-interaction pair(IN0 -> OUT0),
 * every time compute a 1x2 child interaction using AVX
*/
inline void c2c_1x27x1_AVX(real*& M_, real*& IN, real*& OUT)
{
    union{
      __m256d pk = _mm256_set1_pd(0);
      double data[4];
    };
    __m256d qk_riri, qk_irir, gk, gk_0022, gk_1133, pk_a, pk_b;

    qk_riri = _mm256_broadcast_pd((const __m128d*)IN);
    qk_irir = _mm256_permute_pd(qk_riri, 5);

    for(int i = 0; i < 56; i += 4)
    {
        gk = _mm256_load_pd(M_ + i);
        gk_0022 = _mm256_unpacklo_pd(gk, gk);
        gk_1133 = _mm256_unpackhi_pd(gk, gk);
        pk_a = _mm256_mul_pd(qk_riri, gk_0022);
        pk_b = _mm256_mul_pd(qk_irir, gk_1133);
        pk = _mm256_add_pd(pk, _mm256_addsub_pd(pk_a, pk_b));
    }
    pk = _mm256_add_pd(pk, _mm256_permute2f128_pd(pk,pk,1));
    OUT[0] += data[0]; 
    OUT[1] += data[1];
}

/**
 * @brief Compute Pk=Qk*Gk for 8x8 child-interactions in 2 parent-interaction pair(IN0 -> OUT0, IN1 -> OUT1),
 * every time compute two 2x2 child interactions
*/
inline void c2c_8x8x2(real*& M_, real*& IN0, real*& IN1, real*& OUT0, real*& OUT1)
{
    real out_reg000, out_reg001, out_reg010, out_reg011;
    real out_reg100, out_reg101, out_reg110, out_reg111;
    real  in_reg000,  in_reg001,  in_reg010,  in_reg011;
    real  in_reg100,  in_reg101,  in_reg110,  in_reg111;
    real   m_reg000,   m_reg001,   m_reg010,   m_reg011;
    real   m_reg100,   m_reg101,   m_reg110,   m_reg111;
    for(int i1=0;i1<8;i1+=2){
      real* IN0_=IN0;
      real* IN1_=IN1;
      out_reg000=OUT0[ 0]; out_reg001=OUT0[ 1];
      out_reg010=OUT0[ 2]; out_reg011=OUT0[ 3];
      out_reg100=OUT1[ 0]; out_reg101=OUT1[ 1];
      out_reg110=OUT1[ 2]; out_reg111=OUT1[ 3];
      for(int i2=0;i2<8;i2+=2){
        m_reg000=M_[ 0]; m_reg001=M_[ 1];
        m_reg010=M_[ 2]; m_reg011=M_[ 3];
        m_reg100=M_[16]; m_reg101=M_[17];
        m_reg110=M_[18]; m_reg111=M_[19];

        in_reg000=IN0_[0]; in_reg001=IN0_[1];
        in_reg010=IN0_[2]; in_reg011=IN0_[3];
        in_reg100=IN1_[0]; in_reg101=IN1_[1];
        in_reg110=IN1_[2]; in_reg111=IN1_[3];

        out_reg000 += m_reg000*in_reg000 - m_reg001*in_reg001;
        out_reg001 += m_reg000*in_reg001 + m_reg001*in_reg000;
        out_reg010 += m_reg010*in_reg000 - m_reg011*in_reg001;
        out_reg011 += m_reg010*in_reg001 + m_reg011*in_reg000;

        out_reg000 += m_reg100*in_reg010 - m_reg101*in_reg011;
        out_reg001 += m_reg100*in_reg011 + m_reg101*in_reg010;
        out_reg010 += m_reg110*in_reg010 - m_reg111*in_reg011;
        out_reg011 += m_reg110*in_reg011 + m_reg111*in_reg010;

        out_reg100 += m_reg000*in_reg100 - m_reg001*in_reg101;
        out_reg101 += m_reg000*in_reg101 + m_reg001*in_reg100;
        out_reg110 += m_reg010*in_reg100 - m_reg011*in_reg101;
        out_reg111 += m_reg010*in_reg101 + m_reg011*in_reg100;

        out_reg100 += m_reg100*in_reg110 - m_reg101*in_reg111;
        out_reg101 += m_reg100*in_reg111 + m_reg101*in_reg110;
        out_reg110 += m_reg110*in_reg110 - m_reg111*in_reg111;
        out_reg111 += m_reg110*in_reg111 + m_reg111*in_reg110;

        M_+=32; // Jump to (column+2).
        IN0_+=4;
        IN1_+=4;
      }
      OUT0[ 0]=out_reg000; OUT0[ 1]=out_reg001;
      OUT0[ 2]=out_reg010; OUT0[ 3]=out_reg011;
      OUT1[ 0]=out_reg100; OUT1[ 1]=out_reg101;
      OUT1[ 2]=out_reg110; OUT1[ 3]=out_reg111;
      M_+=4-64*2; // Jump back to first column (row+2).
      OUT0+=4;
      OUT1+=4;
    }
  }

/**
 * @brief Compute Pk=Qk*Gk for 8x8 child-interactions in 2 parent-interaction pair(IN0 -> OUT0, IN1 -> OUT1),
 * every time compute two 2x2 child interactions using AVX instructions
*/
inline void c2c_8x8x2_avx(double*& M_, double*& IN0, double*& IN1, double*& OUT0, double*& OUT1)
{
    __m256d out00,out01,out10,out11;
    __m256d out20,out21,out30,out31;
    double* in0__ = IN0;
    double* in1__ = IN1;
    out00 = _mm256_load_pd(OUT0);
    out01 = _mm256_load_pd(OUT1);
    out10 = _mm256_load_pd(OUT0+4);
    out11 = _mm256_load_pd(OUT1+4);
    out20 = _mm256_load_pd(OUT0+8);
    out21 = _mm256_load_pd(OUT1+8);
    out30 = _mm256_load_pd(OUT0+12);
    out31 = _mm256_load_pd(OUT1+12);
    for(int i2=0;i2<8;i2+=2){
      __m256d m00;
      __m256d ot00;
      __m256d mt0,mtt0;
      __m256d in00,in00_r,in01,in01_r;
      in00 = _mm256_broadcast_pd((const __m128d*)in0__);
      in00_r = _mm256_permute_pd(in00,5);
      in01 = _mm256_broadcast_pd((const __m128d*)in1__);
      in01_r = _mm256_permute_pd(in01,5);
      m00 = _mm256_load_pd(M_);
      mt0 = _mm256_unpacklo_pd(m00,m00);
      ot00 = _mm256_mul_pd(mt0,in00);
      mtt0 = _mm256_unpackhi_pd(m00,m00);
      out00 = _mm256_add_pd(out00,_mm256_addsub_pd(ot00,_mm256_mul_pd(mtt0,in00_r)));
      ot00 = _mm256_mul_pd(mt0,in01);
      out01 = _mm256_add_pd(out01,_mm256_addsub_pd(ot00,_mm256_mul_pd(mtt0,in01_r)));
      m00 = _mm256_load_pd(M_+4);
      mt0 = _mm256_unpacklo_pd(m00,m00);
      ot00 = _mm256_mul_pd(mt0,in00);
      mtt0 = _mm256_unpackhi_pd(m00,m00);
      out10 = _mm256_add_pd(out10,_mm256_addsub_pd(ot00,_mm256_mul_pd(mtt0,in00_r)));
      ot00 = _mm256_mul_pd(mt0,in01);
      out11 = _mm256_add_pd(out11,_mm256_addsub_pd(ot00,_mm256_mul_pd(mtt0,in01_r)));
      m00 = _mm256_load_pd(M_+8);
      mt0 = _mm256_unpacklo_pd(m00,m00);
      ot00 = _mm256_mul_pd(mt0,in00);
      mtt0 = _mm256_unpackhi_pd(m00,m00);
      out20 = _mm256_add_pd(out20,_mm256_addsub_pd(ot00,_mm256_mul_pd(mtt0,in00_r)));
      ot00 = _mm256_mul_pd(mt0,in01);
      out21 = _mm256_add_pd(out21,_mm256_addsub_pd(ot00,_mm256_mul_pd(mtt0,in01_r)));
      m00 = _mm256_load_pd(M_+12);
      mt0 = _mm256_unpacklo_pd(m00,m00);
      ot00 = _mm256_mul_pd(mt0,in00);
      mtt0 = _mm256_unpackhi_pd(m00,m00);
      out30 = _mm256_add_pd(out30,_mm256_addsub_pd(ot00,_mm256_mul_pd(mtt0,in00_r)));
      ot00 = _mm256_mul_pd(mt0,in01);
      out31 = _mm256_add_pd(out31,_mm256_addsub_pd(ot00,_mm256_mul_pd(mtt0,in01_r)));
      in00 = _mm256_broadcast_pd((const __m128d*) (in0__+2));
      in00_r = _mm256_permute_pd(in00,5);
      in01 = _mm256_broadcast_pd((const __m128d*) (in1__+2));
      in01_r = _mm256_permute_pd(in01,5);
      m00 = _mm256_load_pd(M_+16);
      mt0 = _mm256_unpacklo_pd(m00,m00);
      ot00 = _mm256_mul_pd(mt0,in00);
      mtt0 = _mm256_unpackhi_pd(m00,m00);
      out00 = _mm256_add_pd(out00,_mm256_addsub_pd(ot00,_mm256_mul_pd(mtt0,in00_r)));
      ot00 = _mm256_mul_pd(mt0,in01);
      out01 = _mm256_add_pd(out01,_mm256_addsub_pd(ot00,_mm256_mul_pd(mtt0,in01_r)));
      m00 = _mm256_load_pd(M_+20);
      mt0 = _mm256_unpacklo_pd(m00,m00);
      ot00 = _mm256_mul_pd(mt0,in00);
      mtt0 = _mm256_unpackhi_pd(m00,m00);
      out10 = _mm256_add_pd(out10,_mm256_addsub_pd(ot00,_mm256_mul_pd(mtt0,in00_r)));
      ot00 = _mm256_mul_pd(mt0,in01);
      out11 = _mm256_add_pd(out11,_mm256_addsub_pd(ot00,_mm256_mul_pd(mtt0,in01_r)));
      m00 = _mm256_load_pd(M_+24);
      mt0 = _mm256_unpacklo_pd(m00,m00);
      ot00 = _mm256_mul_pd(mt0,in00);
      mtt0 = _mm256_unpackhi_pd(m00,m00);
      out20 = _mm256_add_pd(out20,_mm256_addsub_pd(ot00,_mm256_mul_pd(mtt0,in00_r)));
      ot00 = _mm256_mul_pd(mt0,in01);
      out21 = _mm256_add_pd(out21,_mm256_addsub_pd(ot00,_mm256_mul_pd(mtt0,in01_r)));
      m00 = _mm256_load_pd(M_+28);
      mt0 = _mm256_unpacklo_pd(m00,m00);
      ot00 = _mm256_mul_pd(mt0,in00);
      mtt0 = _mm256_unpackhi_pd(m00,m00);
      out30 = _mm256_add_pd(out30,_mm256_addsub_pd(ot00,_mm256_mul_pd(mtt0,in00_r)));
      ot00 = _mm256_mul_pd(mt0,in01);
      out31 = _mm256_add_pd(out31,_mm256_addsub_pd(ot00,_mm256_mul_pd(mtt0,in01_r)));
      M_ += 32;
      in0__ += 4;
      in1__ += 4;
    }
    _mm256_store_pd(OUT0,out00);
    _mm256_store_pd(OUT1,out01);
    _mm256_store_pd(OUT0+4,out10);
    _mm256_store_pd(OUT1+4,out11);
    _mm256_store_pd(OUT0+8,out20);
    _mm256_store_pd(OUT1+8,out21);
    _mm256_store_pd(OUT0+12,out30);
    _mm256_store_pd(OUT1+12,out31);
  }

}