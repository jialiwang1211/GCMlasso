#' Framingham Heart study dataset
#'
#' This dataset contains 4135 patients with 16 risk factors. Some artificial missing
#' data are created for \code{TenYearCHD}.
#'
#' @source Kaggle Framingham Heart study dataset
#' \url{https://www.kaggle.com/amanajmera1/framingham-heart-study-dataset/}
#' @format A data frame with columns:
#' \describe{
#'  \item{male}{binary, male or female.}
#'  \item{age}{continuous, age of the patient.}
#'  \item{education}{categorial, levels coded 1 for some high school, 2 for a high school diploma or GED,
#'   3 for some college or vocational school, and 4 for a college degree.}
#'  \item{currentSmoker}{binary, whether or not the patient is a current smoker.}
#'  \item{cigsPerDay}{continuous, the number of cigarettes that the person smoked on average in one day.}
#'  \item{BPMeds}{binary, whether or not the patient was on blood pressure medication.}
#'  \item{prevalentStroke}{binary, whether or not the patient had previously had a stroke.}
#'  \item{prevalentHyp}{binary, whether or not the patient was hypertensive.}
#'  \item{diabetes}{binary, whether or not the patient had diabetes.}
#'  \item{totChol}{continuous, total cholesterol level.}
#'  \item{sysBP}{continuous, systolic blood pressure.}
#'  \item{diaBP}{continuous, diastolic blood pressure.}
#'  \item{BMI}{continuous, Body Mass Index.}
#'  \item{heartRate}{continuous, heart rate.}
#'  \item{glucose}{continuous, glucose level.}
#'  \item{TenYearCHD}{binary, 10 year risk of coronary heart disease CHD.}
#' }
#' @examples
#' data(Framingham)
#' str(Framingham)
"Framingham"

