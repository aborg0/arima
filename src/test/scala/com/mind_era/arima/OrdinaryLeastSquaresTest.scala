package com.mind_era.arima

import com.mind_era.arima.OrdinaryLeastSquares.{Coefficient, CoefficientErrors, Result, SimpleCoeff, XY}
import org.scalatest._
import org.scalatest.Inside._
import org.scalatest.Matchers._
import org.scalactic._
import spire.math.{Rational, Real}
import spire.implicits._
import scalin.mutable._
import scalin.mutable.dense._
import scalin.syntax._

/**
  * Created by aborg on 04/06/2017.
  */
class OrdinaryLeastSquaresTest extends FreeSpec {

  "Methods tests" - {
    "olsQuadratic" in {

    }

    "ols" in {
      // From https://en.wikipedia.org/w/index.php?title=Ordinary_least_squares&oldid=776500856#Example_with_real_data
      val xs: IndexedSeq[Rational] =
        "1.47\t1.50\t1.52\t1.55\t1.57\t1.60\t1.63\t1.65\t1.68\t1.70\t1.73\t1.75\t1.78\t1.80\t1.83".split('\t').map(
          Rational(_))
      val ys: IndexedSeq[Rational] =
        "52.21\t53.12\t54.48\t55.84\t57.20\t58.57\t59.93\t61.29\t63.11\t64.47\t66.28\t68.10\t69.92\t72.19\t74.46".split(
          '\t').map(Rational(_))
      val result =
        OrdinaryLeastSquares.ols(xs.zip(ys).map((pair) => XY(IndexedSeq(pair._1, pair._1 * pair._1), pair._2)))
      inside(result) {
        case Result(interceptCont, Seq(beta1Cont, beta2Cont), _) =>
          inside(interceptCont) {
            case intercept: Coefficient[_] => assert(intercept.value.toDouble === 128.8128 +- 5e-5)
          }
          inside(beta1Cont) {
            case beta1: Coefficient[_] => assert(beta1.value.toDouble === -143.1620 +- 5e-5)
          }
          inside(beta2Cont) {
            case beta2: Coefficient[_] => assert(beta2.value.toDouble === 61.9603 +- 5e-5)
          }
      }
    }
    "olsWithErrors" in {
      import Conversions._
      val xs: IndexedSeq[Rational] =
        "1.47\t1.50\t1.52\t1.55\t1.57\t1.60\t1.63\t1.65\t1.68\t1.70\t1.73\t1.75\t1.78\t1.80\t1.83".split('\t').map(
          Rational(_))
      val ys: IndexedSeq[Rational] =
        "52.21\t53.12\t54.48\t55.84\t57.20\t58.57\t59.93\t61.29\t63.11\t64.47\t66.28\t68.10\t69.92\t72.19\t74.46".split(
          '\t').map(Rational(_))
      val result: Result[Rational, CoefficientErrors[Rational, Real]] = OrdinaryLeastSquares.olsWithErrors[Rational, Real](xs.zip(ys).map((pair) =>
        XY(IndexedSeq(pair._1, pair._1 * pair._1), pair._2)))
      inside(result) {
        case Result(interceptCont, Seq(beta1Cont, beta2Cont), _) =>
          inside(interceptCont) {
            case CoefficientErrors(intercept, iStdErr, it, ip) => {
              assert(intercept.toDouble === 128.8128 +- 5e-5)
              assert(iStdErr.toDouble === 16.3083 +- 5e-5)
              assert(it.toDouble === 7.8986 +- 5e-5)
              assert(ip.toDouble === 0d +- 5e-5)
            }
          }
          inside(beta1Cont) {
            case CoefficientErrors(beta1, b1StdErr, b1t, b1p) => {
              assert(beta1.toDouble === -143.1620 +- 5e-5)
              assert(b1StdErr.toDouble === 19.8332 +- 5e-5)
              assert(b1t.toDouble === -7.2183 +- 5e-5)
              assert(b1p.toDouble === 0d +- 5e-5)
            }
          }
          inside(beta2Cont) {
            case CoefficientErrors(beta2, b2StdErr, b2t, b2p) => {
              assert(beta2.toDouble === 61.9603 +- 5e-5)
              assert(b2StdErr.toDouble === 6.0084 +- 5e-5)
              assert(b2t.toDouble === 10.3122 +- 5e-5)
              assert(b2p.toDouble === 0d +- 5e-5)
            }
          }
      }
    }

    "olsWithoutInterceptWithErrors" in {
      import Conversions._
      val xs: IndexedSeq[Rational] =
        "1.47\t1.50\t1.52\t1.55\t1.57\t1.60\t1.63\t1.65\t1.68\t1.70\t1.73\t1.75\t1.78\t1.80\t1.83".split('\t').map(
          Rational(_))
      val ys: IndexedSeq[Rational] =
        "52.21\t53.12\t54.48\t55.84\t57.20\t58.57\t59.93\t61.29\t63.11\t64.47\t66.28\t68.10\t69.92\t72.19\t74.46".split(
          '\t').map(Rational(_))
      val result: Result[Rational, CoefficientErrors[Rational, Real]] = OrdinaryLeastSquaresWithoutIntercept.olsWithErrors[Rational, Real](xs.zip(ys).map((pair) =>
        XY(IndexedSeq(pair._1, pair._1 * pair._1), pair._2)))
      inside(result) {
        case Result(interceptCont, Seq(beta1Cont, beta2Cont), _) =>
          inside(interceptCont) {
            case CoefficientErrors(intercept, iStdErr, it, ip) => {
              assert(intercept.toDouble === 0d +- 5e-5)
              assert(iStdErr.toDouble === 0d +- 5e-5)
              assert(it.toDouble === 0d +- 5e-5)
//              assert(ip.toDouble === 0d +- 5e-5)
            }
          }
          inside(beta1Cont) {
            case CoefficientErrors(beta1, b1StdErr, b1t, b1p) => {
              assert(beta1.toDouble === 13.4211 +- 5e-5)
              assert(b1StdErr.toDouble === 1.4323 +- 5e-5)
              assert(b1t.toDouble === 9.37 +- 5e-3)
              assert(b1p.toDouble === 0d +- 5e-5)
            }
          }
          inside(beta2Cont) {
            case CoefficientErrors(beta2, b2StdErr, b2t, b2p) => {
              assert(beta2.toDouble === 14.5867 +- 5e-5)
              assert(b2StdErr.toDouble === 0.8582 +- 5e-5)
              assert(b2t.toDouble === 17d +- 5e-3)
              assert(b2p.toDouble === 0d +- 5e-5)
            }
          }
      }
    }
  }
}
