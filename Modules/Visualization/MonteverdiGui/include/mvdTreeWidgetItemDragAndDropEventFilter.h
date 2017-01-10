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
#ifndef mvdTreeWidgetItemDragAndDropEventFilter_h
#define mvdTreeWidgetItemDragAndDropEventFilter_h

//
// Configuration include.
//// Included at first position before any other ones.
#include "ConfigureMonteverdi.h"


/*****************************************************************************/
/* INCLUDE SECTION                                                           */

//
// Qt includes (sorted by alphabetic order)
//// Must be included before system/custom includes.
#include <QtGui>

//
// System includes (sorted by alphabetic order)

//
// ITK includes (sorted by alphabetic order)

//
// OTB includes (sorted by alphabetic order)
#include "OTBMonteverdiGUIExport.h"
//
// Monteverdi includes (sorted by alphabetic order)
#include "mvdAbstractDragAndDropEventFilter.h"

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


/*****************************************************************************/
/* CLASS DEFINITION SECTION                                                  */

/**
 * \class TreeWidgetItemDragAndDropEventFilter
 *
 * \ingroup OTBMonteverdiGUI
 *
 * \brief Widget template skeleton to copy-paste when adding a new
 * widget class.
 */
class OTBMonteverdiGUI_EXPORT TreeWidgetItemDragAndDropEventFilter :
    public AbstractDragAndDropEventFilter
{

  /*-[ QOBJECT SECTION ]-----------------------------------------------------*/

  Q_OBJECT;

  /*-[ PUBLIC SECTION ]------------------------------------------------------*/

//
// Public methods.
public:

  /** \brief Constructor. */
  TreeWidgetItemDragAndDropEventFilter( QObject* p =NULL );

  /** \brief Destructor. */
  ~TreeWidgetItemDragAndDropEventFilter() ITK_OVERRIDE;

  /*-[ PUBLIC SLOTS SECTION ]------------------------------------------------*/

//
// Public SLOTS.
public slots:

  /*-[ SIGNALS SECTION ]-----------------------------------------------------*/

//
// Signals.
signals:

  /**
   */
  void ItemDropped( QTreeWidgetItem* item );


  /*-[ PROTECTED SECTION ]---------------------------------------------------*/

//
// Protected methods.
protected:

  /**
   * \see http://qt-project.org/doc/qt-4.8/qwidget.html#dragEnterEvent
   */
  bool DragEnterEvent( QObject* object, QDragEnterEvent* event ) ITK_OVERRIDE;

  /**
   * \see http://qt-project.org/doc/qt-4.8/qwidget.html#dragLeaveEvent
   */
  bool DragLeaveEvent( QObject* object, QDragLeaveEvent* event ) ITK_OVERRIDE;

  /**
   * \see http://qt-project.org/doc/qt-4.8/qwidget.html#dragMoveEvent
   */
  bool DragMoveEvent( QObject* object, QDragMoveEvent* event ) ITK_OVERRIDE;

  /**
   * \see http://qt-project.org/doc/qt-4.8/qwidget.html#dropEvent
   */
  bool DropEvent( QObject* object, QDropEvent* event ) ITK_OVERRIDE;

//
// Protected attributes.
protected:

  /*-[ PRIVATE SECTION ]-----------------------------------------------------*/

//
// Private methods.
private:

//
// Private attributes.
private:

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

#endif // mvdTreeWidgetItemDragAndDropEventFilter_h
