/*=========================================================================

  Program:   Monteverdi
  Language:  C++


  Copyright (c) Centre National d'Etudes Spatiales. All rights reserved.
  See Copyright.txt for details.

  Monteverdi is distributed under the CeCILL licence version 2. See
  Licence_CeCILL_V2-en.txt or
  http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt for more details.

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef mvdHistogramPlotPicker_h
#define mvdHistogramPlotPicker_h

//
// Configuration include.
//// Included at first position before any other ones.
#include "ConfigureMonteverdi.h"


/*****************************************************************************/
/* INCLUDE SECTION                                                           */

//
// Qt includes (sorted by alphabetic order)
//// Must be included before system/custom includes.

//
//
#include <qwt_plot_picker.h>

//
// System includes (sorted by alphabetic order)

//
// ITK includes (sorted by alphabetic order)

//
// OTB includes (sorted by alphabetic order)
#include "OTBMonteverdiGUIExport.h"
//
// Monteverdi includes (sorted by alphabetic order)
#include "mvdCore.h"
#include "mvdTypes.h"


/*****************************************************************************/
/* PRE-DECLARATION SECTION                                                   */

//
// External classes pre-declaration.
namespace
{
}

class QwtPlotCanvas;
class QwtPlotCurve;

namespace mvd
{

//
// Internal classes pre-declaration.


/*****************************************************************************/
/* CLASS DEFINITION SECTION                                                  */

/**
 * \class HistogramPlotPicker
 *
 * \ingroup OTBMonteverdiGUI
 *
 * \brief Widget template skeleton to copy-paste when adding a new
 * widget class.
 */
class OTBMonteverdiGUI_EXPORT HistogramPlotPicker :
    public QwtPlotPicker
{

  /*-[ QOBJECT SECTION ]-----------------------------------------------------*/

  Q_OBJECT;

  /*-[ PUBLIC SECTION ]------------------------------------------------------*/

//
// Public types and constants.
public:

  /**
   * \brief
   */
  typedef std::vector< QwtPlotCurve* > PlotCurveVector;

//
// Public methods.
public:

  /** \brief Constructor. */
  HistogramPlotPicker( const PlotCurveVector& curves, QwtPlotCanvas* canvas );

  /** \brief Constructor. */
  HistogramPlotPicker( const PlotCurveVector& curves,
		       int xAxis,
		       int yAxis,
		       QwtPlotCanvas* canvas );

  /** \brief Constructor. */
  HistogramPlotPicker( const PlotCurveVector& curves,
		       int xAxis,
		       int yAxis,
		       DisplayMode trackerMode,
		       QwtPlotCanvas* canvas );

  /** \brief Destructor. */
  ~HistogramPlotPicker() ITK_OVERRIDE;

  /**
   */
  void SetRubberBandPen( RgbwChannel channel, const QPen& pen );

  /**
   */
  void SetGrayscaleActivated( bool activated );

  //
  // QwtPlotPicker methods.

  void drawRubberBand( QPainter* painter ) const ITK_OVERRIDE;

  /*-[ PUBLIC SLOTS SECTION ]------------------------------------------------*/

//
// Public SLOTS.
public slots:

  /*-[ SIGNALS SECTION ]-----------------------------------------------------*/

//
// Signals.
signals:

  /*-[ PROTECTED SECTION ]---------------------------------------------------*/

//
// Protected methods.
protected:

  //
  // QwtPlotPicker methods.

  using QwtPlotPicker::trackerText;

  QwtText trackerText( const QwtDoublePoint & ) const ITK_OVERRIDE;

//
// Protected attributes.
protected:

  /*-[ PRIVATE SECTION ]-----------------------------------------------------*/

//
// Private methods.
private:

  /**
   */
  double Find( const QwtPlotCurve* curve, double x ) const;

  /**
   */
  CountType Find( const QwtPlotCurve* curve,
		  double x,
		  double& xmin,
		  double& xmax,
		  double& y ) const;

//
// Private attributes.
private:
  /**
   * \brief
   */
  static const CountType CURVE_COUNT = 4;

  /**
   * \brief
   */
  PlotCurveVector m_PlotCurves;

  /**
   */
  QPen m_RubberBandPens[ CURVE_COUNT ];

  /**
   */
  bool m_IsGrayscaleActivated : 1;

  /*-[ PRIVATE SLOTS SECTION ]-----------------------------------------------*/

//
// Slots.
private slots:
};

} // end namespace 'mvd'

/*****************************************************************************/
/* INLINE SECTION                                                            */

namespace mvd
{

} // end namespace 'mvd'

#endif // mvdHistogramPlotPicker_h
