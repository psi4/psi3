C     ///////////////////////////////////////////////////////////////
C  TORS     (OPT2)
      SUBROUTINE HIJKS6(NAD,K1,K2,K3,K4,XA,H111,H112,H221,H222,H113,
     $       H123,H223,H331,H332,H333,H411,H421,H422,H431,H432,H433,
     $       H441,H442,H443,H444)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XA(NAD,3),V1(3),V2(3),V3(3),V4(3),E21(3),E23(3),E34(3)
      DIMENSION BP21(3),BP22(3),BP23(3),BP32(3),BP33(3),BP34(3)
      DIMENSION H11(3,3),H21(3,3),H31(3,3),H41(3,3),H22(3,3),H32(3,3)
      DIMENSION H42(3,3),H33(3,3),H43(3,3),H44(3,3)
      DIMENSION H111(3,3,3),H112(3,3,3),H221(3,3,3),H222(3,3,3)
      DIMENSION H113(3,3,3),H123(3,3,3),H223(3,3,3),H331(3,3,3)
      DIMENSION H332(3,3,3),H333(3,3,3),H411(3,3,3),H421(3,3,3)
      DIMENSION H422(3,3,3),H431(3,3,3),H432(3,3,3),H433(3,3,3)
      DIMENSION H441(3,3,3),H442(3,3,3),H443(3,3,3),H444(3,3,3)
      PARAMETER(ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0)
      CALL VECT6(NAD,K1,K2,K3,K4,V1,V2,V3,V4,XA,TAU)
      CALL VECT1(NAD,K1,K2,E21,XA,T21)
      CALL VECT1(NAD,K3,K2,E23,XA,T23)
      CALL VECT1(NAD,K4,K3,E34,XA,T34)
      CALL VECT2(NAD,K1,K2,K3,BP21,BP22,BP23,XA,P2)
      CALL VECT2(NAD,K2,K3,K4,BP32,BP33,BP34,XA,P3)
      CALL MAT1(H32,E23)
      CALL MAT1(H21,E21)
      CALL MAT1(H43,E34)
      C1=ONE/T21
      C2=ONE/T34
      C3=ONE/T23
      C4=DSIN(P2)
      C5=DCOS(P2)
      C6=C5/C4
      C7=DSIN(P3)
      C8=DCOS(P3)
      C9=C8/C7
      C10=ONE/C4**2
      C11=ONE/C7**2
      C12=C1*C1
      C13=C2*C2
      C14=C3*C3
      C15=T21*C3
      C16=T34*C3
      W1=TWO*C6
      W2=TWO*C9
      W3=TWO*C1
      W4=TWO*C2
      W5=C5*C3
      W6=C8*C3
      DO 10  K=1,3
      H411(1,1,K)=E21(K)*C1+BP21(K)*C6
      H411(1,2,K)=E34(K)*C2+BP34(K)*C9
      H411(1,3,K)=E23(K)*C3+BP23(K)*W1
      H411(2,1,K)=-E23(K)*C3+BP32(K)*W2
      H411(2,2,K)=E21(K)*W3+E23(K)*C3-BP22(K)*W1
      H411(2,3,K)=E34(K)*W4-E23(K)*C3-BP33(K)*W2
      H411(3,1,K)=E23(K)*W5+BP23(K)*C4
 10   H411(3,2,K)=-E23(K)*W6+BP32(K)*C7
      DO 15  K=1,3
        DO 16  M=1,3
 16     V2(M)=ZERO
      V2(K)=ONE
      CALL MAT1(H22,V2)
      DO 20  J=1,3
      DO 20  I=1,3
 20   H421(I,J,K)=H22(I,J)
 15   CONTINUE
      W1=TWO*C10*C12
      W2=TWO*C11*C13
      DO 100 K=1,3
      W3=W1*H411(1,1,K)
      W4=W2*H411(1,2,K)
        DO 100  J=1,K
        DO 100  I=1,J
        H111(I,J,K)=W3*H32(I,J)
 100    H444(I,J,K)=-W4*H32(I,J)
      W1=C10*C12
      W2=C11*C13
      W3=W1*C3
      W4=W2*C3
      DO 105 K=1,3
      W5=W1*H411(1,3,K)
      W6=W2*H411(2,1,K)
        DO 105  J=1,3
        DO 105  I=1,3
        H113(I,J,K)=W5*H32(I,J)-W3*H421(I,J,K)
 105    H442(I,J,K)=-W6*H32(I,J)-W4*H421(I,J,K)
      W1=C1*C3*C10
      W2=C2*C3*C11
      DO 110 K=1,3
        W5=W1*H411(2,2,K)
        W6=W2*H411(2,3,K)
        DO 110  J=1,3
        DO 110  I=1,3
        H123(I,K,J)=-W5*H21(I,J)+W3*H421(I,J,K)
 110    H432(I,K,J)=-W6*H43(I,J)+W4*H421(I,J,K)
      CALL HIJS2(NAD,K1,K2,K3,XA,H11,H21,H31,H22,H32,H33)
      CALL HIJS1(NAD,K1,K2,XA,H44)
      CALL HIJS1(NAD,K2,K3,XA,H42)
      W1=TWO*C1
      W2=TWO*C12
      DO 115  K=1,3
        DO 116  I=1,K
 116    H43(I,K)=TWO*(W1*H44(I,K)+C6*H11(I,K)-C10*BP21(I)*BP21(K))
 115  H43(K,K)=H43(K,K)-W2
      DO 120  K=1,3
      DO 120  J=1,K
      DO 120  I=1,J
 120  H111(I,J,K)=H111(I,J,K)-V1(J)*H43(I,K)
        W1=TWO*C6
        W2=TWO*C10
        W3=TWO*C3
        DO 125  K=1,3
        DO 125  I=1,3
 125    H43(I,K)=H31(K,I)*W1-BP21(I)*BP23(K)*W2
        DO 130  K=1,3
        DO 130  J=1,3
        DO 130  I=1,3
 130    H113(I,J,K)=H113(I,J,K)-V1(J)*H43(I,K)
        DO 135  K=1,3
        DO 136  J=1,3
 136    H43(J,K)=W3*H42(J,K)-W1*H32(K,J)+W2*BP22(J)*BP23(K)
 135    H43(K,K)=H43(K,K)-C14
        DO 140  K=1,3
        DO 140  J=1,3
        DO 140  I=1,3
 140    H123(I,J,K)=H123(I,J,K)+V1(I)*H43(J,K)
      W1=C4*C3
      W2=C4*C15
      W3=C5*C15
      W4=W3*C3
      W5=C3*C15
      DO 145  K=1,3
      DO 145  I=1,3
      W6=-E21(I)*BP23(K)*W1+H32(K,I)*W2+BP22(I)*BP23(K)*W3-H42(I,K)*W4
 145  H43(I,K)=W6
      DO 150  K=1,3
      DO 150  I=1,3
      W6=H43(I,K)+E23(I)*W5*H411(3,1,K)
      DO 150  J=1,3
 150  H223(I,J,K)=-V1(J)*W6
      W1=C3*C4*C15
      W2=C5*C14
      DO 155  K=1,3
      DO 155  I=1,3
      W3=-E23(K)*BP22(I)*W1+E23(K)*W2*(C15*E23(I)-E21(I))
 155  H43(I,K)=H43(I,K)+W3
      DO 160  K=1,3
      DO 160  J=1,3
      DO 160  I=1,3
 160  H332(I,J,K)=V1(J)*H43(K,I)
      CALL HIJS2(NAD,K2,K3,K4,XA,H22,H32,H42,H33,H43,H44)
      CALL HIJS1(NAD,K4,K3,XA,H11)
      CALL HIJS1(NAD,K3,K2,XA,H31)
      W1=TWO*C2
      W2=TWO*C13
      DO 165  K=1,3
      DO 166  I=1,K
 166  H21(I,K)=TWO*(W1*H11(I,K)+C9*H44(I,K)-C11*BP34(I)*BP34(K))
 165  H21(K,K)=H21(K,K)-W2
      DO 170  K=1,3
      DO 170  J=1,K
      DO 170  I=1,J
 170  H444(I,J,K)=H444(I,J,K)-V4(J)*H21(I,K)
      W1=TWO*C9
      W2=TWO*C11
      W3=TWO*C3
      DO 175  K=1,3
      DO 175  I=1,3
 175  H21(I,K)=W1*H42(I,K)-W2*BP34(I)*BP32(K)
      DO 180  K=1,3
      DO 180  J=1,3
      DO 180  I=1,3
 180  H442(I,J,K)=H442(I,J,K)-V4(J)*H21(I,K)
      DO 185  K=1,3
      DO 186  J=1,3
 186  H21(J,K)=W3*H31(J,K)-W1*H32(J,K)+W2*BP33(J)*BP32(K)
 185  H21(K,K)=H21(K,K)-C14
      DO 190  K=1,3
      DO 190  J=1,3
      DO 190  I=1,3
 190  H432(I,J,K)=H432(I,J,K)+V4(I)*H21(J,K)
      W1=C7*C3
      W2=C7*C16
      W3=C8*C16
      W4=W3*C3
      W5=T34*C14
      DO 195  K=1,3
      DO 195  I=1,3
      W6=-E34(I)*BP32(K)*W1+H32(I,K)*W2+BP33(I)*BP32(K)*W3-H31(I,K)*W4
 195  H21(I,K)=W6
      DO 200  K=1,3
      DO 200  I=1,3
      W6=H21(I,K)-E23(I)*W5*H411(3,2,K)
      DO 200  J=1,3
 200  H332(I,J,K)=H332(I,J,K)-V4(J)*W6
      W1=C3*C7*C16
      W2=C8*C14
      DO 205  K=1,3
      DO 205  I=1,3
      W3=E23(K)*BP33(I)*W1+E23(K)*W2*(E34(I)+C16*E23(I))
 205  H21(I,K)=H21(I,K)+W3
      DO 210  K=1,3
      DO 210  J=1,3
      DO 210  I=1,3
 210  H223(I,J,K)=H223(I,J,K)+V4(J)*H21(K,I)
      CALL HIJS6(NAD,K1,K2,K3,K4,XA,H11,H21,H31,H41,H22,H32,H42,
     $            H33,H43,H44)
      DO 220  K=1,3
      DO 220  J=1,K
      DO 220  I=1,J
      H111(I,J,K)=H111(I,J,K)-TWO*H11(J,K)*H411(1,1,I)
 220  H444(I,J,K)=H444(I,J,K)-TWO*H44(J,K)*H411(1,2,I)
      CALL FILL3A(3,3,H111)
      CALL FILL3A(3,3,H444)
      DO 225  I=1,3
      W1=TWO*H411(1,1,I)
      W2=TWO*H411(1,2,I)
      DO 225  J=1,3
      DO 225  K=1,3
      H113(I,J,K)=H113(I,J,K)-W1*H31(K,J)
      H442(I,J,K)=H442(I,J,K)-W2*H42(J,K)
      H123(I,J,K)=H123(I,J,K)-H21(J,I)*H411(1,3,K)
 225  H432(I,J,K)=H432(I,J,K)-H43(I,J)*H411(2,1,K)
      W4=C5*C15
      W1=W4-ONE
      W2=C8*C16
      W3=W2-ONE
      DO 230  K=1,3
      DO 230  J=1,3
      DO 230  I=1,3
      H223(I,J,K)=H223(I,J,K)+W1*H123(J,I,K)-W2*H432(J,K,I)
 230  H332(I,J,K)=H332(I,J,K)+W3*H432(J,I,K)-W4*H123(J,K,I)
      DO 240  K=1,3
      DO 240  J=1,3
      DO 240  I=1,3
      H223(I,J,K)=H223(I,J,K)-C15*H21(I,J)*H411(3,1,K)
 240  H332(I,J,K)=H332(I,J,K)-C16*H43(J,I)*H411(3,2,K)
      DO 250  K=1,3
      DO 250  J=1,3
      W1=C16*(H43(J,K)-C3*V4(J)*E23(K))
      W2=C15*(H21(K,J)+C3*V1(J)*E23(K))
      DO 250  I=1,3
      H223(I,J,K)=H223(I,J,K)+W1*H411(3,2,I)
 250  H332(I,J,K)=H332(I,J,K)+W2*H411(3,1,I)
      W1=C5*C3
      W2=C4*C15
      W3=C5*T21*C14
      W4=C8*C3
      W5=C7*C16
      W6=C8*T34*C14
      DO 260  K=1,3
      H411(1,1,K)=W5*BP33(K)+W6*E23(K)+W4*E34(K)
      H411(1,2,K)=W2*BP22(K)-W3*E23(K)+W1*E21(K)
      H411(1,3,K)=-W1*E21(K)-W2*BP22(K)+W3*E23(K)
 260  H411(2,1,K)=-W4*E34(K)-W5*BP33(K)-W6*E23(K)
      DO 265  K=1,3
      DO 265  J=1,3
      DO 265  I=1,3
      H223(I,J,K)=H223(I,J,K)+H42(J,I)*H411(1,1,K)
 265  H332(I,J,K)=H332(I,J,K)+H31(I,J)*H411(1,2,K)
      DO 270  K=1,3
      DO 270  J=1,3
      W1=(H31(K,J)-C3*V1(J)*E23(K))
      W2=(H42(J,K)+C3*V4(J)*E23(K))
      DO 270  I=1,3
      H223(I,J,K)=H223(I,J,K)+W1*H411(1,3,I)
 270  H332(I,J,K)=H332(I,J,K)+W2*H411(2,1,I)
      DO 280  K=1,3
      DO 280  J=1,3
      DO 280  I=1,3
      H411(I,J,K)=ZERO
      H421(I,J,K)=ZERO
      H431(I,J,K)=ZERO
      H441(I,J,K)=ZERO
      H443(J,K,I)=-(H444(I,J,K)+H442(J,K,I))
 280  H112(I,J,K)=-(H111(I,J,K)+H113(I,J,K))
      DO 285  K=1,3
      DO 285  J=1,3
      DO 285  I=1,3
      H433(K,I,J)=-(H443(I,K,J)+H432(K,J,I))
 285  H221(I,J,K)=-(H112(I,K,J)+H123(K,J,I))
      DO 290  K=1,3
      DO 290  J=1,3
      DO 290  I=1,3
      H422(K,I,J)=-(H432(K,I,J)+H442(I,K,J))
 290  H331(I,J,K)=-(H123(K,I,J)+H113(I,K,J))
      DO 300  K=1,3
      DO 300  J=1,K
      DO 300  I=1,J
      H222(I,J,K)=-(H221(I,J,K)+H223(I,J,K)+H422(K,I,J))
 300  H333(I,J,K)=-(H331(I,J,K)+H332(I,J,K)+H433(K,I,J))
      CALL FILL3A(3,3,H222)
      CALL FILL3A(3,3,H333)
      RETURN
      END
