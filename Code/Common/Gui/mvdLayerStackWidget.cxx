/*=========================================================================

  Program:   Monteverdi2
  Language:  C++


  Copyright (c) Centre National d'Etudes Spatiales. All rights reserved.
  See Copyright.txt for details.

  Monteverdi2 is distributed under the CeCILL licence version 2. See
  Licence_CeCILL_V2-en.txt or
  http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt for more details.

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#include "mvdLayerStackWidget.h"
#include "ui_mvdLayerStackWidget.h"


/*****************************************************************************/
/* INCLUDE SECTION                                                           */

//
// Qt includes (sorted by alphabetic order)
//// Must be included before system/custom includes.

//
// System includes (sorted by alphabetic order)
#include <cassert>

//
// ITK includes (sorted by alphabetic order)

//
// OTB includes (sorted by alphabetic order)

//
// Monteverdi includes (sorted by alphabetic order)
#include "Gui/mvdLayerStackItemModel.h"

namespace mvd
{

/*
  TRANSLATOR mvd::LayerStackWidget

  Necessary for lupdate to be aware of C++ namespaces.

  Context comment for translator.
*/


/*****************************************************************************/
/* CONSTANTS                                                                 */


/*****************************************************************************/
/* STATIC IMPLEMENTATION SECTION                                             */


/*****************************************************************************/
/* CLASS IMPLEMENTATION SECTION                                              */

/*******************************************************************************/
LayerStackWidget
::LayerStackWidget( QWidget* parent, Qt::WindowFlags flags  ):
  QWidget( parent, flags ),
  m_UI( new mvd::Ui::LayerStackWidget() )
{
  m_UI->setupUi( this );

  {
  QItemSelectionModel * ism = m_UI->treeView->selectionModel();

  m_UI->treeView->setModel( new LayerStackItemModel( m_UI->treeView ) );

  delete ism;
  ism = NULL;
  }
}

/*******************************************************************************/
LayerStackWidget
::~LayerStackWidget()
{
  delete m_UI;
  m_UI = NULL;
}

/*******************************************************************************/
const LayerStackItemModel *
LayerStackWidget
::GetItemModel() const
{
  return const_cast< LayerStackWidget * >( this )->GetItemModel();
}

/*******************************************************************************/
LayerStackItemModel *
LayerStackWidget
::GetItemModel()
{
  assert(
    m_UI->treeView->model()==
    qobject_cast< LayerStackItemModel * >( m_UI->treeView->model() )
    );

  return qobject_cast< LayerStackItemModel * >( m_UI->treeView->model() );
}

/*******************************************************************************/
// void
// LayerStackWidget
// ::SetModel( LayerStackItemModel * itemModel )
// {
//   // See http://qt-project.org/doc/qt-4.8/qabstractitemview.html#setModel .
//   QItemSelectionModel * itemSelectionModel = m_UI->treeView->selectionModel();

//   m_UI->treeView->setModel( itemModel );

//   itemModel->setParent( m_UI->treeView );

//   delete itemSelectionModel;
//   itemSelectionModel = NULL;
// }

/*******************************************************************************/
void
LayerStackWidget
::SetSelection( int row )
{
  assert( m_UI->treeView->selectionModel()!=NULL );

  m_UI->treeView->selectionModel()->select(
    m_UI->treeView->model()->index( row, 1 ),
    QItemSelectionModel::ClearAndSelect | QItemSelectionModel::Rows
  );
}

/*******************************************************************************/
/* SLOTS                                                                       */
/*******************************************************************************/

} // end namespace 'mvd'
