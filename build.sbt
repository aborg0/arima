name := "arima"

version := "0.0.1"

scalaVersion := "2.12.6"

moduleName := "com.github.aborg0"

resolvers ++= Seq(
  "bintray/denisrosset/maven" at "https://dl.bintray.com/denisrosset/maven",
  Resolver.sonatypeRepo("snapshots")
)

libraryDependencies += "net.alasc" %% "scalin-core" % "0.16.0.0"

libraryDependencies += "org.scalactic" %% "scalactic" % "3.0.4"

libraryDependencies += "org.scalatest" %% "scalatest" % "3.0.4" % "test"

//libraryDependencies += "org.scalaz" %% "scalaz-core" % "7.2.12"

libraryDependencies += "org.typelevel" %% "cats" % "0.9.0"

libraryDependencies += "com.outr" %% "scribe" % "1.4.3"

scalacOptions ++= Seq("-unchecked", "-deprecation", "-feature")