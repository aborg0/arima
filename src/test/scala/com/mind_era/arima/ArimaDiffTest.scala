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
      val pi = Vec[Rational](3, 1, 4, 1, 5, 9)
      assert (Vec[Rational](-2, 3, -3, 4, 4) === Arima.diff(pi, Natural.one, Natural.one))
    }

    "diff matrix" in {

    }

  }
}
