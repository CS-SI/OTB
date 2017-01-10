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
#ifndef mvdColorSetupController_h
#define mvdColorSetupController_h

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
// System includes (sorted by alphabetic order)

//
// ITK includes (sorted by alphabetic order)

//
// OTB includes (sorted by alphabetic order)
#include "OTBMonteverdiGUIExport.h"
//
// Monteverdi includes (sorted by alphabetic order)
#include "mvdCore.h"
//
#include "mvdAbstractModelController.h"


/*****************************************************************************/
/* PRE-DECLARATION SECTION                                                   */

//
// External classes pre-declaration.
namespace
{
}

namespace mvd
{
//
// Internal classes pre-declaration.
class ColorSetupWidget;

/*****************************************************************************/
/* CLASS DEFINITION SECTION                                                  */

/**
 * \class ColorSetupController
 *
 * \ingroup OTBMonteverdiGUI
 *
 * \brief Color-setup widget controller for VectorImageModel objects.
 */
class OTBMonteverdiGUI_EXPORT ColorSetupController :
    public AbstractModelController
{

  /*-[ QOBJECT SECTION ]-----------------------------------------------------*/

  Q_OBJECT;

  /*-[ PUBLIC SECTION ]------------------------------------------------------*/

//
// Public methods.
public:

  /**
   * \brief Constructor.
   *
   * \param widget Controlled widget.
   * \param parent Parent QObject of this QObject.
   */
  ColorSetupController( ColorSetupWidget* widget, QObject* p =NULL );

  /**
   * \brief Destructor.
   */
  ~ColorSetupController() ITK_OVERRIDE;


  /*-[ SIGNALS SECTION ]-----------------------------------------------------*/

//
// Signals.
signals:
  /**
   * \brief Signal emitted when the band-index of a RGB channel has
   * been changed. This signal is used to forward events to other
   * controllers such as ColorDynamicsController.
   *
   * \param channel The RGB channel for which the band-index has changed.
   * \param band The new band-index of the RGB channel.
   */
  void RgbChannelIndexChanged( RgbwChannel channel, int band );

  /**
   * \brief Signal emitted when the band-index of a white (gray)
   * channel has been changed. This signal is used to forward events
   * to other controllers such as ColorDynamicsController.
   *
   * \param band The new band-index of the white (gray) channel.
   */
  void GrayChannelIndexChanged( int band );

  /**
   * \brief Signal emitted when the grayscale-mode activation-state
   * has changed. This signal is used to forward events to other
   * controllers such as ColorDynamicsController.
   *
   * \param activated The new grayscale-mode activation-state.
   */
  void GrayscaleActivated( bool activated );

  /*-[ PROTECTED SECTION ]---------------------------------------------------*/

//
// Protected methods.
protected:

//
// Protected attributes.
protected:

  /*-[ PRIVATE SECTION ]-----------------------------------------------------*/

//
// Private methods.
private:

  /**
   * \brief Reset current band-index to default values for given RGB
   * channels.
   *
   * \param channels Given channels for which to reset current-band
   * index. \see RgbBound() for valid values.
   */
  void ResetIndices( RgbwChannel channels );

  /**
   */
  void ResetAlpha();

  //
  // AbstractModelController methods.

  void Connect( AbstractModel* ) ITK_OVERRIDE;

  void ClearWidget() ITK_OVERRIDE;

  void virtual_ResetWidget( bool ) ITK_OVERRIDE;

  void Disconnect( AbstractModel* ) ITK_OVERRIDE;

//
// Private attributes.
private:

  /*-[ PRIVATE SLOTS SECTION ]-----------------------------------------------*/

//
// Slots.
private slots:
  /**
   * \brief Slot called when the current index of a RGB channel of the
   * controlled widget has changed.
   *
   * \param channel The RGB video-channel for which the current
   * band-index has changed.
   *
   * \param index The new band-index of the RGB video-channel.
   */
  void OnCurrentRgbIndexChanged( RgbwChannel channel, int index );

  /**
   * \brief Slot called when the current index of the white (gray)
   * channel of the controlled widget has changed.
   *
   * \param index The new band-index of the white (gray)
   * video-channel.
   */
  void OnCurrentGrayIndexChanged( int index );

  /**
   * \brief Slot called when the grayscale-mode has been
   * activated/deactivated.
   *
   * \param activated The new activation state.
   */
  void OnGrayscaleActivated( bool activated );

  /**
   */
  void OnAlphaValueChanged( double );

};

} // end namespace 'mvd'.

/*****************************************************************************/
/* INLINE SECTION                                                            */

namespace mvd
{
} // end namespace 'mvd'

#endif // mvdColorSetupController_h
