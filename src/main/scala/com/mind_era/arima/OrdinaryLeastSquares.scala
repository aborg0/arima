package com.mind_era.arima

import algebra.ring.{Field, Ring}

import scala.math.Numeric._
import scalin.{DenseMat, Pivot}
import scalin.algos.Inverse
import scalin.mutable._
import spire.algebra.{Eq, NRoot}
import spire.implicits._

import scala.math.ScalaNumber
//import spire.syntax.ring._

/**
  * Created by aborg on 04/06/2017.
  *
  * Based on https://en.wikipedia.org/w/index.php?title=Ordinary_least_squares&oldid=776500856
  */
object OrdinaryLeastSquares {
  case class XY[V: Eq: Field](x: IndexedSeq[V], y: V)
  object XY {
    def apply[V: Eq: Field](x: V, y: V): XY[V] = {
      XY(SingleValue(x), y)
    }
  }

  sealed trait Coefficient[V] {
    def value: V
  }

  case class Result[V: Field, C <: Coefficient[V]](intercept: C,
                                                       beta: IndexedSeq[C]/*d*/,
                                                       diff/*epsilon*/: IndexedSeq[V]/*n*/)

  final case class SimpleCoeff[V: Field](value: V) extends Coefficient[V]

  final case class CoefficientErrors[V: Field, ErrV: Field: NRoot](value: V, stdError: ErrV, t: ErrV, p: Double)
    extends Coefficient[V]

  final case class NDimCoefficient[V: Eq: Field: Coefficient](coefficients: IndexedSeq[V])

  private def createMatX[V: Field](values: IndexedSeq[XY[V]]): scalin.mutable.DenseMat[V] = DenseMat.tabulate(
    values.length, 1 + values(0).x.length)((row, col) => if (col == 0) Ring.one else values(row).x(col - 1))

  def ols[V: Field: MatEngine: Pivot: VecEngine](values: IndexedSeq[XY[V]] /*length n*/)(
    implicit denseInverse: Inverse[V, Mat[V]]): Result[V, SimpleCoeff[V]] = {
    val BetaDiff(beta, diff, _) = olsBetaDiff(values)

    Result(
      SimpleCoeff(beta(0, 0)),
      beta = beta(::, 0).toIndexedSeq.tail.map(SimpleCoeff(_)),
      diff = diff(::, 0).toIndexedSeq
    )
  }

  case class BetaDiff[V](beta: Mat[V], diff: Mat[V], `(xᵀx)⁻¹`: Mat[V])

  private def olsBetaDiff[V: Field : MatEngine : Pivot : VecEngine](values: IndexedSeq[XY[V]])(
    implicit denseInverse: Inverse[V, Mat[V]]): BetaDiff[V] = {
    val matX: Mat[V] = createMatX(values)
    val matXT: Mat[V] = matX.t
    val y: DenseVec[V] = DenseVec(values.map(_.y): _*)
    val xᵀx: Mat[V] = matXT * matX
    val `(xᵀx)⁻¹`: Mat[V] = xᵀx.inverse
    val tmp: Mat[V] = `(xᵀx)⁻¹` * matXT
    val beta: Mat[V] = tmp * y.toColMat
    val yHat = matX * beta
    val diff = y.toColMat - yHat
    BetaDiff(beta, diff, `(xᵀx)⁻¹`)
  }

  def olsWithErrors[V: Field: MatEngine: Pivot: VecEngine, ErrV <: ScalaNumber: Field: NRoot](values: IndexedSeq[XY[V]])(
    implicit inverse: Inverse[V, Mat[V]], conv: Convert[V, ErrV]): Result[V, CoefficientErrors[V, ErrV]] = {
    val BetaDiff(beta, diff, xTxInverse) = olsBetaDiff(values)
    val `(xᵀx)⁻¹`: Mat[V] = xTxInverse
    val diffVec: Vec[V] = diff(::, 0)
    val nV: V = Ring.fromInt[V](values.length - values(0).x.length - 1)
    // val sSquare = diff.t * diff / (values.length - values(0).x.length)
    val sigmaHatSquare: V = diffVec.dot(diffVec) / nV
    def sigmaHatJ(j: Int): ErrV = {
      val xTxJ: V = xTxInverse(j, j)
      val square: V = Field.times(sigmaHatSquare, xTxJ)
      conv.to(square).sqrt
    }

    val interceptT = conv.to(beta(0, 0)) / sigmaHatJ(0)
    Result(CoefficientErrors[V, ErrV](beta(0, 0), sigmaHatJ(0), interceptT, pValue(interceptT, sigmaHatJ(0))),
      beta = beta(::, 0).toIndexedSeq.zipWithIndex.tail.map(tuple => {
        val i: Int = tuple._2
        val t = conv.to(beta(i, 0)) / sigmaHatJ(i)
        CoefficientErrors[V, ErrV](tuple._1, sigmaHatJ(i), t, pValue(t, sigmaHatJ(i)))
      }
    ), diff = diff(::, 0).toIndexedSeq)
  }

  def olsQuadratic[V: Eq: Field, Res<: Coefficient[V]](values: IndexedSeq[XY[V]]): Result[V, Res] = {
    ???
  }

  def phiApprox(x: Double): Double = {
      .5 * (1 + math.signum(x)*math.sqrt(- math.expm1(-2/math.Pi * x * x)))
  }

  def erfApprox(x: Double): Double = {
    phiApprox(math.sqrt(2) * x) * 2 - 1
  }

  def pValue[ErrV: Numeric](t: ErrV, sigma: ErrV): Double = {
    val num = implicitly[Numeric[ErrV]]
    erfApprox(num.mkNumericOps(t).toDouble() / 2 * num.mkNumericOps(sigma).toDouble())
  }

  def pValue[ErrV <: ScalaNumber](t: ErrV, sigma: ErrV): Double = {
    erfApprox(math.abs(t.doubleValue()) / 2 * sigma.doubleValue()) - 1d
  }
}


