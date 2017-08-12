name := "arima"

version := "0.0.1"

scalaVersion := "2.12.3"

moduleName := "com.github.aborg0"

libraryDependencies += "net.alasc" %% "scalin-core" % "0.14.1.1-SNAPSHOT"

libraryDependencies += "org.scalactic" %% "scalactic" % "3.0.1"

libraryDependencies += "org.scalatest" %% "scalatest" % "3.0.1" % "test"

//libraryDependencies += "org.scalaz" %% "scalaz-core" % "7.2.12"

libraryDependencies += "org.typelevel" %% "cats" % "0.9.0"

libraryDependencies += "com.outr" %% "scribe" % "1.4.2"

scalacOptions ++= Seq("-unchecked", "-deprecation", "-feature")