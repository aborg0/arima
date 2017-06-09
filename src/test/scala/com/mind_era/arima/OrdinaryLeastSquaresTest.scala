package com.mind_era.arima

import com.mind_era.arima.OrdinaryLeastSquares.{Result, SimpleCoeff, XY}
import org.scalatest.FreeSpec
import org.scalatest.Matchers._
import org.scalactic._
import spire.math.Rational
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
      val Result(SimpleCoeff(intercept), Seq(SimpleCoeff(beta1), SimpleCoeff(beta2)), _) =
        OrdinaryLeastSquares.ols(xs.zip(ys).map((pair) => XY(IndexedSeq(pair._1, pair._1 * pair._1), pair._2)))
      assert(intercept.toDouble === 128.8128 +- 5e-5)
      assert(beta1.toDouble === -143.1620 +- 5e-5)
      assert(beta2.toDouble === 61.9603 +- 5e-5)
    }

  }
}
