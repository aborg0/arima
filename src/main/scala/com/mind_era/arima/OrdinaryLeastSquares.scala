package com.mind_era.arima

import algebra.ring.{Field, Ring}
import scalin.{DenseMat, Pivot, mutable}
import scalin.algos.Inverse
import scalin.mutable._
import spire.algebra.Eq

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

  case class Result[V: Eq: Field, C <: Coefficient[V]](intercept: C,
                                                       beta: IndexedSeq[C]/*d*/,
                                                       diff/*epsilon*/: IndexedSeq[V]/*n*/)

  final case class SimpleCoeff[V: Eq: Field](value: V) extends Coefficient[V]

  final case class CoefficientErrors[V: Eq: Field](value: V, stdError: V, t: V, p: V) extends Coefficient[V]

  final case class NDimCoefficient[V: Eq: Field: Coefficient](coefficients: IndexedSeq[V])

  private def createMatX[V: Eq: Field](values: IndexedSeq[XY[V]]): scalin.mutable.DenseMat[V] = scalin.mutable.DenseMat.tabulate(values.length,
    1 + values(0).x.length)((row, col) => if (col == 0) Ring.one else values(row).x(col - 1))

  def ols[V: Eq: Field: MatEngine: Pivot: VecEngine](values: IndexedSeq[XY[V]] /*length n*/)(implicit denseInverse: Inverse[V, Mat[V]]): Result[V, SimpleCoeff[V]] = {
    //implicit val denseMatEngine: scalin.mutable.MatEngine[V] = scalin.mutable.DenseMat.defaultEngine
    val matX: scalin.mutable.Mat[V] = createMatX(values)
    val matXT: scalin.mutable.Mat[V] = matX.t[scalin.mutable.Mat[V]](scalin.mutable.DenseMat.defaultEngine)
    val y: scalin.mutable.DenseVec[V] = DenseVec(values.map(_.y): _*)
//    implicit val denseMatEngine: scalin.mutable.MatEngine[V] = scalin.mutable.DenseMat.defaultEngine
    val xᵀx: scalin.mutable.Mat[V] = matXT * matX
    val `(xᵀx)⁻¹`: scalin.mutable.Mat[V] = xᵀx.inverse
    val tmp: scalin.mutable.Mat[V] = `(xᵀx)⁻¹` * matXT
    val beta: Mat[V] = tmp * y.toColMat
    val yHat = matX * beta
    val diff = y.toColMat - yHat
    Result(SimpleCoeff(beta(0, 0)), beta = beta(::, 0).toIndexedSeq.tail.map(SimpleCoeff(_)), diff = diff(::, 0).toIndexedSeq)
  }

  def olsQuadratic[V: Eq: Field, Res<: Coefficient[V]](values: IndexedSeq[XY[V]]): Result[V, Res] = {
    ???
  }
}


