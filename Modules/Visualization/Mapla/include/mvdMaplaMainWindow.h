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

#ifndef mvdMaplaMainWindow_h
#define mvdMaplaMainWindow_h

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

//
// Monteverdi includes (sorted by alphabetic order)
#include "mvdTypes.h"
//
#include "mvdI18nMainWindow.h"

#include "OTBMaplaExport.h"

/*****************************************************************************/
/* PRE-DECLARATION SECTION                                                   */

//
// External class pre-declaration.
namespace
{
}

namespace mvd
{
//
// Internal class pre-declaration.
// Core
// Gui
class ApplicationsToolBoxController;

namespace Ui
{
class MaplaMainWindow;
}

/*****************************************************************************/
/* CLASS DEFINITION SECTION                                                  */

/**
 * \class MaplaMainWindow
 *
 * \ingroup OTBMapla
 *
 * \brief The application main-widow widget.
 */
class OTBMapla_EXPORT MaplaMainWindow
  : public I18nMainWindow
{

  /*-[ QOBJECT SECTION ]-----------------------------------------------------*/

  Q_OBJECT;

  /*-[ PUBLIC SECTION ]------------------------------------------------------*/

//
// Public types.
public:

//
// Public methods.
public:
  /** \brief Constructor. */
  MaplaMainWindow( QWidget* Parent =0, Qt::WindowFlags flags =0 );

  /** \brief Destructor. */
  ~MaplaMainWindow() ITK_OVERRIDE;

  /*-[ SIGNALS SECTION ]-----------------------------------------------------*/

//
// SIGNALS.
signals:

  /*-[ PROTECTED SECTION ]---------------------------------------------------*/

//
// Protected methods.
protected:

  //
  // QMainWindow methods.

  void closeEvent( QCloseEvent* event ) ITK_OVERRIDE;

//
// Protected attributes.
protected:

  /*-[ PROTECTED SLOTS SECTION ]---------------------------------------------*/

//
// Protected slots.
protected slots:

  /**
   */
  void OnAboutToChangeModel( const AbstractModel* ) ITK_OVERRIDE;

  /**
   */
  void OnModelChanged( AbstractModel* ) ITK_OVERRIDE;

  /**
   */
  void OnApplicationToLaunchSelected( const QString & name,
                                      const QString & doc );

  /**
   */
#if 0
  void OnOTBApplicationOutputImageChanged( const QString & name,
                                           const QString & filename);
#endif
  /**
   */
  void OnExecutionDone( int status );

  /*-[ PRIVATE SECTION ]-----------------------------------------------------*/

//
// Private types.
private:


//
// Private methods.
private:

  /**
   */
  void InitializeDockWidgets();

  /**
   */
  void InitializeCentralWidget();

  //
  // I18nMainWindow methods.

  void virtual_SetupUI() ITK_OVERRIDE;

  void virtual_ConnectUI() ITK_OVERRIDE;

//
// Private attributes.
private:
  /**
   * \brief uic generated.
   */
  Ui::MaplaMainWindow* m_UI;

#ifdef OTB_USE_QT4
  /**
   *  \brief OTB-applications tool-box controller.
   */
  ApplicationsToolBoxController * m_ApplicationsToolBoxController;
#endif

  /*-[ PRIVATE SLOTS SECTION ]-----------------------------------------------*/

//
// Private slots.
private slots:
  /**
   * \brief Qt auto-connected slot which is called when
   * 'File/Preferences' menu action is trigerred.
   */
  void on_action_Preferences_triggered();
};

} // end namespace 'mvd'

/*****************************************************************************/
/* INLINE SECTION                                                            */

//
// System includes (sorted by alphabetic order)

//
// ITK includes (sorted by alphabetic order)

//
// OTB includes (sorted by alphabetic order)

//
// Monteverdi includes (sorted by alphabetic order)
#include "mvdImageViewWidget.h"

//
// Some constants.

namespace mvd
{

} // end namespace 'mvd'

#endif // maplaMainWindow_h
