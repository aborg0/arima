package com.mind_era.arima

import org.scalatest.FreeSpec
//import scalin.immutable.Vec
//import scalin._
import spire.algebra._
//import cats.kernel.Eq._
import spire.math._
import spire.implicits._
import scalin.mutable._
import scalin.mutable.dense._

class ArimaDiffTest extends FreeSpec {

  "Methods tests" - {
    "diff vector" in {
//      implicit val eq = Eq.fromUniversalEquals[Number]
//      implicit val field: Field[Int] = Field.
      // > s <- c(3, 1, 4, 1, 5, 9)
      val pi = Vec[Rational](3, 1, 4, 1, 5, 9)
      // > diff(s, 1, 1)
      // [1] -2  3 -3  4  4
      assert (Vec[Rational](-2, 3, -3, 4, 4) === Arima.diff(pi, Natural.one, Natural.one))
      // > diff(s, 1, 3)
      // [1] -11  13  -7
      assert (Vec[Rational](-11, 13, -7) === Arima.diff(pi, Natural.one, Natural(3)))
      // > diff(s, 2, 1)
      // [1] 1 0 1 8
      assert (Vec[Rational](1, 0, 1, 8) === Arima.diff(pi, Natural(2), Natural.one))
    }

    "diff matrix" in {
      // > m <- matrix(s, 2, 3)
      //      [,1] [,2] [,3]
      // [1,]    3    4    5
      // [2,]    1    1    9
      val pi23 = Mat[Rational]((3: Rational, 4: Rational, 5: Rational), (1: Rational, 1: Rational, 9: Rational))
      // > diff(m, 1, 1)
      //      [,1] [,2] [,3]
      // [1,]   -2   -3    4
      assert (Vec[Rational](-2, -3, 4).toRowMat === Arima.diff(pi23, Natural.one, Natural.one))
      // > diff(m, 1, 2)
      // numeric(0)
      assert (Mat.fillConstant(0, 3)(Rational.zero) === Arima.diff(pi23, Natural.one, Natural(2)))
      // > diff(m, 2, 1)
      // numeric(0)
      assert (Mat.fillConstant(0, 3)(Rational.zero) === Arima.diff(pi23, Natural(2), Natural.one))
      // 3.14159265359
      // > mm <- matrix(c(3,1,4,1,5,9,2,6,5,3,5,9), 3, 4)
      //      [,1] [,2] [,3] [,4]
      // [1,]    3    1    2    3
      // [2,]    1    5    6    5
      // [3,]    4    9    5    9
      val row1: (Rational, Rational, Rational, Rational) = (3, 1, 2, 3)
      val row2: (Rational, Rational, Rational, Rational) = (1, 5, 6, 5)
      val row3: (Rational, Rational, Rational, Rational) = (4, 9, 5, 9)
      val pi34 = Mat(row1, row2, row3)
      // > diff(mm, 2, 1)
      //      [,1] [,2] [,3] [,4]
      // [1,]    1    8    3    6
      assert (Vec[Rational](1, 8, 3, 6).toRowMat === Arima.diff(pi34, Natural(2), Natural.one))
      // > diff(mm, 1, 2)
      //      [,1] [,2] [,3] [,4]
      // [1,]    5    0   -5    2
      assert (Vec[Rational](5, 0, -5, 2).toRowMat === Arima.diff(pi34, Natural.one, Natural(2)))

      // > diff(mm, 1, 1)
      //      [,1] [,2] [,3] [,4]
      // [1,]   -2    4    4    2
      // [2,]    3    4   -1    4
      val resRow1: (Rational, Rational, Rational, Rational) = (-2, 4, 4, 2)
      val resRow2: (Rational, Rational, Rational, Rational) = (3, 4, -1, 4)
      assert (Mat(resRow1, resRow2) === Arima.diff(pi34, Natural.one, Natural.one))
    }

  }
}
