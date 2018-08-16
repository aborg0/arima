package com.mind_era.arima

import com.mind_era.arima.Arima.SVD
import scalin.Mat
import scalin.Vec
import spire.algebra.{Eq, Field}
import spire.std.DoubleAlgebra
import scalin.Mat.eye
import scalin.mutable.DenseMat
import scalin.mutable.dense._
import scalin.syntax.assign._
import spire.implicits._ // for the spire.std instances

class ComputeSVD[V: Eq: Field : DoubleAlgebra]() {
  // TODO implement
  def svd(matrix: Mat[V]): SVD[V] = ???
  case class Output1a(B: Mat[V], U: Mat[V], V: Mat[V])
  def algorithm1a(m: Int, n: Int, matrix: Mat[V]): Output1a = {
    val U: Mat[V] = irregularEye(m, n)
    val V: Mat[V] = eye[V](n)
    (0 until n).foreach{k =>
      val Qk = getQk(matrix.colSeq[V].apply(k).to[Vec], k - 1)
      matrix := Qk * matrix
      U := U * Qk
    }
  }

  def irregularEye (m: Int, n: Int): DenseMat[V] = {
    val zeroMatrix = Mat.fill[V](m, n)(Field.zero[V])
    val identityMatrix = eye[V](Math.min(m, n))
    zeroMatrix + identityMatrix
  }

  def getQk(column: Vec[V], k: Int): Mat[V] = {
    val s: V = column dot column
    val v: Vec[V] = Vec.fill[V](column.length)(Field.zero)
    val m: Vec[V] = Vec.fill[V](column.length)(Field.zero)
    v(k - 1) := column(k-1)
    v(k) := s
    (k - 1 until column.length).foreach(i => m(i) := column(i))
    ???
  }

}
