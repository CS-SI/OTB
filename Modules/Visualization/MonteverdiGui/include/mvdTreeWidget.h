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
#ifndef mvdTreeWidget_h
#define mvdTreeWidget_h

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


#define USE_CUSTOM_MIME_DATA 1


/*****************************************************************************/
/* PRE-DECLARATION SECTION                                                   */

//
// External classes pre-declaration.
namespace
{
}

namespace mvd
{

/*****************************************************************************/
/* FUNCTIONS DEFINITION SECTION                                              */

/**
 */
QMimeData*
EncodeMimeData( QMimeData* mimeData, const QList< QTreeWidgetItem* >& items );

/**
 */
int
DecodeMimeData( QList< QTreeWidgetItem* >& items, const QMimeData* mimeData );

//
// Internal classes pre-declaration.
namespace Ui
{
//class TreeWidget;
};


/*****************************************************************************/
/* CLASS DEFINITION SECTION                                                  */

/**
 * \class TreeWidget
 *
 * \ingroup OTBMonteverdiGUI
 *
 * \brief Widget template skeleton to copy-paste when adding a new
 * widget class.
 */
class OTBMonteverdiGUI_EXPORT TreeWidget :
    public QTreeWidget
{

  /*-[ QOBJECT SECTION ]-----------------------------------------------------*/

  Q_OBJECT;

  /*-[ PUBLIC SECTION ]------------------------------------------------------*/

//
// Public types and constants.
public:
  /**
   */
  static const char* ITEM_MIME_TYPE;

  /**
   */
  typedef QList< QTreeWidgetItem* > QTreeWidgetItemList;

//
// Public methods.
public:

  /** \brief Constructor. */
  TreeWidget( QWidget* p =NULL );

  /** \brief Destructor. */
  ~TreeWidget() ITK_OVERRIDE;

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
  void ItemMoved( QTreeWidgetItem * item, QTreeWidgetItem * target );

  /*-[ PROTECTED SECTION ]---------------------------------------------------*/

//
// Protected methods.
protected:

  //
  // QTreeWidget overloads.

  QStringList mimeTypes() const ITK_OVERRIDE;

  QMimeData* mimeData( const QList< QTreeWidgetItem* > items ) const ITK_OVERRIDE;

  void dragEnterEvent( QDragEnterEvent* event ) ITK_OVERRIDE;
  void dragMoveEvent( QDragMoveEvent* event ) ITK_OVERRIDE;
  void dragLeaveEvent( QDragLeaveEvent* event ) ITK_OVERRIDE;
  void dropEvent( QDropEvent* event ) ITK_OVERRIDE;

  Qt::DropActions supportedDropActions() const ITK_OVERRIDE;
  void startDrag( Qt::DropActions supportedActions ) ITK_OVERRIDE;

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
/* GLOBAL FUNCTIONS DECLARATION                                              */

#if USE_CUSTOM_MIME_DATA

//
// Declare Qt tree-widget item pointer types so they can be wrapped by
// QVariant.
Q_DECLARE_METATYPE( QTreeWidgetItem* );


#define TREE_WIDGET_ITEM_USE_STREAM_OPERATORS 1

#if TREE_WIDGET_ITEM_USE_STREAM_OPERATORS

/**
 */
QDataStream&
operator << ( QDataStream& out, QTreeWidgetItem const * item );

/**
 */
QDataStream&
operator >>( QDataStream& in, QTreeWidgetItem * & item );

#endif // !DATASTREAM_USE_TEMPLATE_OPERATORS

#endif // USE_CUSTOM_MIME_DATA

/*****************************************************************************/
/* INLINE SECTION                                                            */

namespace mvd
{
} // end namespace 'mvd'

#endif // mvdTreeWidget_h
