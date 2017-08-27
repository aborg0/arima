package com.mind_era.arima

import org.scalactic.TypeCheckedTripleEquals._
import org.scalatest.{FreeSpec, Inside}
import scalin._
import scalin.mutable.dense._
import spire.implicits._
import spire.math.{Natural, Real}
import Natural._
import com.mind_era.arima.Arima.PhiTheta
import org.scalactic.{Equality, TolerantNumerics}
import scalin.mutable.Vec

class ArimaTransformTest extends FreeSpec with Inside {
//> library(stats)
//> debug(arima)
//> arima(c(3.1, 4.1, 5.9, 2.1, 0, 0), c(1, 0, 1))
  implicit val doubleEquality: Equality[Double] = TolerantNumerics.tolerantDoubleEquality(1E-5)

  "Methods tests" - {
    "invParTrans" - {
      "p == 1" in {
        val coeffs: IndexedSeq[Real] = IndexedSeq[Real](0, 0, 59)
        //Browse[2]> .Call(C_ARIMA_Invtrans, as.double(c(0, 0, 59)), as.integer(c(2, 0, 0,0,0,0)))
        //[1]  0  0 59
        assert(coeffs === Arima.invParTrans(one, coeffs))
        //Browse[2]> .Call(C_ARIMA_Invtrans, as.double(c(0.5, 0, .59)), as.integer(c(1, 0, 0,0,0,0)))
        //  [1] 0.5493061 0.0000000 0.5900000
        inside(Arima.invParTrans(one, coeffs.updated(0,.5: Real))) {
          case Seq(first, second, third) =>
            assert(0.5493061 === first.toDouble)
            assert(0d === second.toDouble)
            assert(59 === third.toDouble)
        }
      }
      "p == 2" in {
        //Browse[2]> .Call(C_ARIMA_Invtrans, as.double(c(0.5, 0.25, .9)), as.integer(c(2, 0, 0,0,0,0)))
        //  [1] 0.8047190 0.2554128 0.9000000
        val coeffs: IndexedSeq[Real] = IndexedSeq[Real](0.5, 0.25, .9)
        inside(Arima.invParTrans(Natural(2), coeffs)) {
          case Seq(first, second, third) =>
            assert(0.8047190 === first.toDouble)
            assert(0.2554128 === second.toDouble)
            assert(.9 === third.toDouble)
        }
      }
    }

    "undoPars" in {
      //Browse[2]> .Call(C_ARIMA_undoPars,
      // as.double(c( 0.8047190, 0.2554128, 0.2209164, 0.2766926, 0.1256572, 0.125, 0.25)),
      // as.integer(c(2, 0, 3, 0, 0, 0)))
      //  [1] 0.500 0.250 0.125 0.250 0.125 0.125 0.250
      val coeffs: IndexedSeq[Real] = IndexedSeq[Real](0.5, 0.25, .125, .25, .125, .125, .25)
      val input = Arima.arimaInvTrans(coeffs, IndexedSeq(Natural(2), zero, Natural(3), zero))
      inside(Arima.undoPars(input.toIndexedSeq
        /*IndexedSeq(
          0.8047190: Real, 0.2554128: Real, 0.2209164: Real, 0.2766926: Real, 0.1256572: Real,
          0.125: Real, 0.25: Real)*/, IndexedSeq(Natural(2), zero, Natural(3), zero, zero, zero)).toIndexedSeq) {
        case IndexedSeq(p0, p1, p2, p3, p4, p5, p6) =>
          assert (p0.toDouble === .5)
          assert (p1.toDouble === .25)
          assert (p2.toDouble === .125)
          assert (p3.toDouble === .25)
          assert (p4.toDouble === .125)
          assert (p5.toDouble === .125)
          assert (p6.toDouble === .25)
      }

    }

    "transPars" in {
      //Browse[2]> .Call(C_ARIMA_transPars,
      // as.double(c( 0.8047190, 0.2554128, 0.2209164, 0.2766926, 0.1256572, 0.125, 0.25)),
      // as.integer(c(2, 0, 0, 0, 0, 0)), TRUE)
      //  [[1]]
      //[1] 0.50 0.25
      //
      //  [[2]]
      //numeric(0)
      val coeffs: IndexedSeq[Real] = IndexedSeq[Real](0.5, 0.25, .125, .25, .125, .125, .25)
      val arma: IndexedSeq[Natural] = IndexedSeq(Natural(2), zero, Natural(3), zero, zero, zero)
      val input: IndexedSeq[Real] = Arima.arimaInvTrans(coeffs, arma).toIndexedSeq
      inside (Arima.transPars(input, arma, true)) {
        case PhiTheta(phi, theta) =>
          assert(phi(0) === .5)
          assert(phi(1) === .25)
          assert(theta === IndexedSeq.empty)
      }
    }

    "arimaInvTrans" - {// See also invParTrans
      "leading zeroes" in {
        val p = Natural(2)
        val coeffs: IndexedSeq[Real] = IndexedSeq[Real](0, 0, .59)
        assert(Vec(coeffs: _*) === Arima.arimaInvTrans(coeffs, IndexedSeq(zero, zero, zero, zero)))
        assert(Vec(coeffs: _*) === Arima.arimaInvTrans(coeffs, IndexedSeq(zero, one, zero, zero)))
        // atanh of 0 is 0, so the results are the same.
        assert(Vec(coeffs: _*) === Arima.arimaInvTrans(coeffs, IndexedSeq(one, one, zero, zero)))
        assert(Vec(coeffs: _*) === Arima.arimaInvTrans(coeffs, IndexedSeq(p, zero, zero, zero)))
        assert(Vec(coeffs: _*) === Arima.arimaInvTrans(coeffs, IndexedSeq(p, one, zero, zero)))
      }
      "also seasonal p" in {
        val p = Natural(2)
        val coeffs: IndexedSeq[Real] = IndexedSeq[Real](0.5, 0.25, .125, .25, .125, .125, .25)
        //Browse[2]> .Call(C_ARIMA_Invtrans, as.double(c(0.5, 0.25, .125, .25, .125, .125, .25)),
        //                 as.integer(c(1, 0, 3, 0, 0, 0)))
        //  [1] 0.5493061 0.3942287 0.2027326 0.2554128 0.1250000 0.1250000 0.2500000
        inside (Arima.arimaInvTrans(coeffs, IndexedSeq(one, zero, Natural(3), zero))) {
          case vec =>
            assert(coeffs.drop(4) === vec.toIndexedSeq.drop(4))
            assert(vec(0).toDouble === 0.5493061)
            assert(vec(1).toDouble === 0.3942287)
            assert(vec(2).toDouble === 0.2027326)
            assert(vec(3).toDouble === 0.2554128)
        }
        //Browse[2]> .Call(C_ARIMA_Invtrans, as.double(c(0.5, 0.25, .125, .25, .125, .125, .25)),
        // as.integer(c(0, 0, 3, 0, 0, 0)))
        //  [1] 1.0732904 0.3288202 0.1256572 0.2500000 0.1250000 0.1250000 0.2500000
        inside (Arima.arimaInvTrans(coeffs, IndexedSeq(zero, zero, Natural(3), zero))) {
          case vec =>
            assert(coeffs.drop(3) === vec.toIndexedSeq.drop(3))
            assert(vec(0).toDouble === 1.0732904)
            assert(vec(1).toDouble === 0.3288202)
            assert(vec(2).toDouble === 0.1256572)
        }
        //Browse[2]> .Call(C_ARIMA_Invtrans, as.double(c(0.5, 0.25, .125, .25, .125, .125, .25)),
        // as.integer(c(2, 0, 3, 0, 0, 0)))
        //  [1] 0.8047190 0.2554128 0.2209164 0.2766926 0.1256572 0.1250000 0.2500000
        inside (Arima.arimaInvTrans(coeffs, IndexedSeq(Natural(2), zero, Natural(3), zero))) {
          case vec =>
            assert(coeffs.drop(5) === vec.toIndexedSeq.drop(5))
            assert(vec(0).toDouble === 0.8047190)
            assert(vec(1).toDouble === 0.2554128)
            assert(vec(2).toDouble === 0.2209164)
            assert(vec(3).toDouble === 0.2766926)
            assert(vec(4).toDouble === 0.1256572)
        }
      }
    }

    "parTrans" - {
      "p == 1" in {
        //Browse[2]> .Call(C_ARIMA_undoPars,
        // as.double(c( 0.8047190, 0.2554128, 0.2209164, 0.2766926, 0.1256572, 0.125, 0.25)),
        // as.integer(c(1, 0, 0, 0, 0, 0)))
        //  [1] 0.6666667 0.2554128 0.2209164 0.2766926 0.1256572 0.1250000 0.2500000
        inside(Arima.parTrans(one,
          IndexedSeq(
            0.8047190: Real, 0.2554128: Real, 0.2209164: Real, 0.2766926: Real, 0.1256572: Real,
            0.125: Real, 0.25: Real))) {
          case IndexedSeq(p0, p1, p2, p3, p4, p5, p6) =>
            assert (p0.toDouble === .66666666)
            assert (p1.toDouble === 0.2554128)
            assert (p2.toDouble === 0.2209164)
            assert (p3.toDouble === 0.2766926)
            assert (p4.toDouble === 0.1256572)
            assert (p5.toDouble === .125)
            assert (p6.toDouble === .25)
        }
      }
      "p == 2" in {
        //Browse[2]> .Call(C_ARIMA_undoPars,
        // as.double(c( 0.8047190, 0.2554128, 0.2209164, 0.2766926, 0.1256572, 0.125, 0.25)),
        // as.integer(c(2, 0, 0, 0, 0, 0)))
        //  [1] 0.5000000 0.2500000 0.2209164 0.2766926 0.1256572 0.1250000 0.2500000
        inside(Arima.parTrans(Natural(2),
          IndexedSeq(
            0.8047190: Real, 0.2554128: Real, 0.2209164: Real, 0.2766926: Real, 0.1256572: Real,
            0.125: Real, 0.25: Real))) {
          case IndexedSeq(p0, p1, p2, p3, p4, p5, p6) =>
            assert (p0.toDouble === .5)
            assert (p1.toDouble === .25)
            assert (p2.toDouble === 0.2209164)
            assert (p3.toDouble === 0.2766926)
            assert (p4.toDouble === 0.1256572)
            assert (p5.toDouble === .125)
            assert (p6.toDouble === .25)
        }
      }
    }

  }
}
